#include "edge_table.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "edge_format.hpp"
#include "edge_shuffle.hpp"
#include "edge_text.hpp"
#include "py_repr.hpp"

using namespace edge_text;

namespace {

// How pandas' C parser would read one value of a column it infers.
enum class Kind { kNa, kInt, kFloat, kStr, kUnsure };

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// [+-]?[0-9]+ ; `digits` = significant digits
bool int_token(std::string_view t, int &digits) {
    size_t p = (!t.empty() && (t[0] == '+' || t[0] == '-')) ? 1 : 0;
    if (p == t.size()) return false;
    digits = 0;
    bool lead = true;
    for (; p < t.size(); ++p) {
        if (!is_digit(t[p])) return false;
        if (lead && t[p] == '0') continue;
        lead = false;
        ++digits;
    }
    return true;
}

// [+-]? (D+ (. D*)? | . D+) ([eE] [+-]? D+)? ; `digits` = significant
// mantissa digits
bool float_token(std::string_view t, int &digits) {
    size_t p = (!t.empty() && (t[0] == '+' || t[0] == '-')) ? 1 : 0;
    digits = 0;
    int mant = 0;
    bool lead = true;
    auto take = [&] {
        for (; p < t.size() && is_digit(t[p]); ++p, ++mant) {
            if (lead && t[p] == '0') continue;
            lead = false;
            ++digits;
        }
    };
    take();
    if (p < t.size() && t[p] == '.') {
        ++p;
        take();
    }
    if (mant == 0) return false;
    if (p < t.size() && (t[p] == 'e' || t[p] == 'E')) {
        ++p;
        if (p < t.size() && (t[p] == '+' || t[p] == '-')) ++p;
        const size_t e0 = p;
        while (p < t.size() && is_digit(t[p])) ++p;
        if (p == e0) return false;
    }
    return p == t.size();
}

Kind classify(std::string_view t) {
    if (is_na_token(t)) return Kind::kNa;
    const char c0 = t[0];
    int digits = 0;
    if (is_digit(c0) || c0 == '+' || c0 == '-' || c0 == '.') {
        if (int_token(t, digits)) return digits <= 18 ? Kind::kInt : Kind::kUnsure;
        if (float_token(t, digits)) return digits <= 17 ? Kind::kFloat : Kind::kUnsure;
        return Kind::kUnsure;  // e.g. "-inf", "1e", "1,5"
    }
    if (std::isspace(static_cast<unsigned char>(c0)) ||
        std::isspace(static_cast<unsigned char>(t.back())))
        return Kind::kUnsure;
    std::string low(t);
    for (auto &c : low) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (low == "inf" || low == "infinity" || low == "nan" || low == "true" || low == "false")
        return Kind::kUnsure;  // float or bool to pandas
    return Kind::kStr;
}

struct ColumnInfo {
    bool na = false, i = false, f = false, s = false;
    // the dtype pandas infers: object, float64 (also int with NA), int64
    Kind kind() const {
        if (s) return Kind::kStr;
        if (f || na) return Kind::kFloat;
        if (i) return Kind::kInt;
        return Kind::kFloat;  // no values at all: NaN
    }
};

// Value fields (columns 2..) of a stored row, at most n_values of them;
// returns how many the row has (0 when it had only the gene columns).
int split_rest(std::string_view rest, int n_values, std::string_view *out) {
    if (rest.empty() || n_values <= 0) return 0;
    const char *p = rest.data(), *e = p + rest.size();
    int n = 0;
    while (n < n_values) {
        const char *t = static_cast<const char *>(std::memchr(p, '\t', static_cast<size_t>(e - p)));
        const char *fe = t ? t : e;
        out[n++] = std::string_view(p, static_cast<size_t>(fe - p));
        if (!t) break;
        p = t + 1;
    }
    return n;
}

double to_double(std::string_view t) {
    std::string tmp(t);
    return std::strtod(tmp.c_str(), nullptr);
}

void write_all(std::FILE *fh, std::string &buf, const std::string &path) {
    if (!buf.empty() && std::fwrite(buf.data(), 1, buf.size(), fh) != buf.size()) {
        std::fclose(fh);
        throw std::runtime_error("error writing " + path);
    }
    buf.clear();
}

// The table pandas would read: genes interned, the value fields of each row
// kept as text, and the dtype pandas infers for each value column.
struct EdgeTable {
    std::vector<std::string> names;
    std::vector<int64_t> a, b;   // gene codes per row
    std::string rest;            // value fields of every row, back to back
    std::vector<uint64_t> rest_off{0};
    std::vector<Kind> kinds;     // per value column
    std::string why;             // non-empty: not handled

    std::string_view row_rest(int64_t k) const {
        return std::string_view(rest.data() + rest_off[k],
                                static_cast<size_t>(rest_off[k + 1] - rest_off[k]));
    }
};

// Reads spec.path; with `keep`, only the rows whose first or second gene
// (unquoted) is in it are kept. The column dtypes and the checks come from
// every row (read_network reads the whole file), except that a missing gene
// name only matters in a kept row.
EdgeTable read_edge_table(const EdgeTableSpec &spec,
                          const std::unordered_set<std::string_view> *keep) {
    EdgeTable t;
    const int ncol = spec.n_columns;
    if (ncol < 2) throw std::invalid_argument("a network needs at least two columns");
    const int nval = ncol - 2;
    std::unordered_map<std::string, int64_t> ids;
    std::vector<ColumnInfo> cols(static_cast<size_t>(nval));
    std::vector<const char *> fs(static_cast<size_t>(ncol) + 1), fe(fs.size());
    std::string sa, sb, sv;
    bool header_pending = spec.skip_header;
    int64_t line_no = 0;

    auto id_of = [&](std::string_view s) {
        auto it = ids.find(std::string(s));
        if (it != ids.end()) return it->second;
        const auto id = static_cast<int64_t>(t.names.size());
        t.names.emplace_back(s);
        ids.emplace(t.names.back(), id);
        return id;
    };
    auto fail = [&](const std::string &what) {
        t.why = "line " + std::to_string(line_no) + ": " + what;
    };

    for_each_line(spec.path, [&](const char *lb, const char *le) {
        ++line_no;
        if (!t.why.empty()) return;
        if (header_pending) { header_pending = false; return; }
        if (lb == le) return;  // pandas skips blank lines
        const int n = split_tabs(lb, le, fs.data(), fe.data(), ncol + 1);
        const std::string_view ga = unquote(fs[0], fe[0], sa);
        const std::string_view gb = n > 1 ? unquote(fs[1], fe[1], sb) : std::string_view();
        if (n > ncol) return fail("more fields than columns");
        const bool kept = !keep || keep->count(ga) || (n > 1 && keep->count(gb));
        if (kept && (n < 2 || is_na_token(ga) || is_na_token(gb)))
            return fail("missing gene name");
        for (int c = 0; c < nval; ++c) {  // dtypes: every row, as pandas reads them all
            ColumnInfo &ci = cols[static_cast<size_t>(c)];
            if (2 + c >= n) { ci.na = true; continue; }
            switch (classify(unquote(fs[2 + c], fe[2 + c], sv))) {
                case Kind::kNa: ci.na = true; break;
                case Kind::kInt: ci.i = true; break;
                case Kind::kFloat: ci.f = true; break;
                case Kind::kStr: ci.s = true; break;
                case Kind::kUnsure: return fail("value pandas may read differently");
            }
        }
        if (!kept) return;
        t.a.push_back(id_of(ga));
        t.b.push_back(id_of(gb));
        if (n > 2) t.rest.append(fs[2], le);
        t.rest_off.push_back(t.rest.size());
    });
    if (!t.why.empty()) return t;
    for (const auto &ci : cols) t.kinds.push_back(ci.kind());
    if (spec.v1_layout) {
        if (nval != 2 && nval != 3) throw std::invalid_argument("v1 layout needs 4 or 5 columns");
        for (const Kind k : t.kinds)
            if (k == Kind::kStr) t.why = "non-numeric correlation, p-value or FDR";
    }
    if (t.names.size() > (size_t(1) << 31)) throw std::invalid_argument("too many genes");
    return t;
}

// Writes spec.header and the rows of `t`, with gene_b codes `b`, as
// Network.to_edgelist does.
void write_edge_table(const EdgeTableSpec &spec, const EdgeTable &t,
                      const std::vector<int64_t> &b) {
    const int nval = spec.n_columns - 2;
    const int64_t n_edges = static_cast<int64_t>(t.a.size());
    std::FILE *fh = std::fopen(spec.output.c_str(), "wb");
    if (!fh) throw std::runtime_error("cannot write " + spec.output);
    std::string buf = spec.header;
    constexpr int64_t kChunk = 1 << 20;
    if (spec.v1_layout) {
        // names as is, correlation %f, p-value / FDR %e
        write_all(fh, buf, spec.output);
        const bool with_q = nval == 3;
        for (int64_t k0 = 0; k0 < n_edges; k0 += kChunk) {
            const int64_t k1 = std::min(n_edges, k0 + kChunk);
            buf = format_edge_lines(
                t.names, static_cast<size_t>(k1 - k0), with_q, spec.n_threads, [&](size_t i) {
                    const int64_t k = k0 + static_cast<int64_t>(i);
                    std::string_view v[3];
                    const int n = split_rest(t.row_rest(k), nval, v);
                    thread_local std::string scratch;
                    double x[3];
                    for (int c = 0; c < 3; ++c) {
                        x[c] = std::numeric_limits<double>::quiet_NaN();
                        if (c < n) {
                            const std::string_view u =
                                unquote(v[c].data(), v[c].data() + v[c].size(), scratch);
                            if (!is_na_token(u)) x[c] = to_double(u);
                        }
                    }
                    return EdgeRow{static_cast<int>(t.a[k]), static_cast<int>(b[k]), x[0], x[1],
                                   x[2]};
                });
            write_all(fh, buf, spec.output);
        }
    } else {
        // DataFrame.to_csv(sep="\t", index=False) of the inferred dtypes
        std::vector<std::string_view> v(static_cast<size_t>(std::max(nval, 1)));
        std::string scratch;
        for (int64_t k = 0; k < n_edges; ++k) {
            append_csv_field(buf, t.names[t.a[k]]);
            buf += '\t';
            append_csv_field(buf, t.names[b[k]]);
            const int n = split_rest(t.row_rest(k), nval, v.data());
            for (int c = 0; c < nval; ++c) {
                buf += '\t';
                if (c >= n) continue;  // missing: NaN, written empty
                const std::string_view u = unquote(v[c].data(), v[c].data() + v[c].size(), scratch);
                if (is_na_token(u)) continue;
                switch (t.kinds[static_cast<size_t>(c)]) {
                    case Kind::kInt: {
                        std::string tmp(u);
                        buf += std::to_string(std::strtoll(tmp.c_str(), nullptr, 10));
                        break;
                    }
                    case Kind::kFloat: append_py_repr(buf, to_double(u)); break;
                    default: append_csv_field(buf, u); break;
                }
            }
            buf += '\n';
            if (buf.size() >= (size_t(1) << 22)) write_all(fh, buf, spec.output);
        }
        write_all(fh, buf, spec.output);
    }
    if (std::fclose(fh) != 0) throw std::runtime_error("error writing " + spec.output);
}

}  // namespace

EdgeTableResult shuffle_edge_file(const EdgeTableSpec &spec, int64_t swap_factor, Pcg64 &rng) {
    EdgeTableResult res;
    EdgeTable t = read_edge_table(spec, nullptr);
    if (!t.why.empty()) {
        res.reason = t.why;
        return res;
    }
    const int64_t n_edges = static_cast<int64_t>(t.a.size());
    res.edge_count = n_edges;
    res.node_count = static_cast<int64_t>(t.names.size());
    std::vector<int64_t> b;  // gene_b codes after the swaps (gene_a never changes)
    {
        EdgeSwapper sw(t.a, std::move(t.b), std::max<int64_t>(1, res.node_count));
        if (n_edges >= 2) res.swaps_applied = sw.shuffle(rng, swap_factor * n_edges);
        b = sw.b();
    }
    write_edge_table(spec, t, b);
    res.ok = true;
    return res;
}

EdgeTableResult extract_edge_file(const EdgeTableSpec &spec,
                                  const std::vector<std::string> &genes) {
    EdgeTableResult res;
    const std::unordered_set<std::string_view> keep(genes.begin(), genes.end());
    const EdgeTable t = read_edge_table(spec, &keep);
    if (!t.why.empty()) {
        res.reason = t.why;
        return res;
    }
    res.edge_count = static_cast<int64_t>(t.a.size());
    res.node_count = static_cast<int64_t>(t.names.size());
    write_edge_table(spec, t, t.b);
    res.ok = true;
    return res;
}
