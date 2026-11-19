#include "edge_files.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "edge_text.hpp"
#include "parallel.hpp"
#include "py_repr.hpp"

using namespace edge_text;

MergedEdges merge_edge_files(const std::vector<EdgeFileSpec> &files, double cor_cutoff) {
    std::unordered_map<std::string, int32_t> ids;
    std::vector<std::string> names;
    auto id_of = [&](std::string_view s) -> int32_t {
        auto it = ids.find(std::string(s));
        if (it != ids.end()) return it->second;
        const auto id = static_cast<int32_t>(names.size());
        names.emplace_back(s);
        ids.emplace(names.back(), id);
        return id;
    };
    // Per ordered pair: Kahan sum as in pandas' group_sum.
    std::unordered_map<uint64_t, int64_t> pair_index;
    std::vector<uint64_t> pair_key;
    std::vector<double> sum, comp;
    std::string sa, sb, sc;

    for (const auto &f : files) {
        if (f.cor_column < 2) throw std::invalid_argument("cor_column must be >= 2");
        const int want = f.cor_column + 1;
        std::vector<const char *> fs(want), fe(want);
        bool header_pending = f.skip_header;
        int64_t line_no = 0;
        for_each_line(f.path, [&](const char *b, const char *e) {
            ++line_no;
            if (header_pending) { header_pending = false; return; }
            if (b == e) return;  // blank lines are skipped
            const int n = split_tabs(b, e, fs.data(), fe.data(), want);
            if (n < 2) return;  // gene_b missing: pandas drops the NaN group key
            const std::string_view a = unquote(fs[0], fe[0], sa);
            const std::string_view g = unquote(fs[1], fe[1], sb);
            if (is_na_token(a) || is_na_token(g)) return;
            const double v = n > f.cor_column
                                 ? parse_value(unquote(fs[f.cor_column], fe[f.cor_column], sc),
                                               f.path, line_no)
                                 : std::numeric_limits<double>::quiet_NaN();
            const int32_t ia = id_of(a);
            const int32_t ib = id_of(g);
            const uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(ia)) << 32) |
                                 static_cast<uint32_t>(ib);
            auto ins = pair_index.emplace(key, static_cast<int64_t>(pair_key.size()));
            if (ins.second) {
                pair_key.push_back(key);
                sum.push_back(0.0);
                comp.push_back(0.0);
            }
            if (std::isnan(v)) return;  // skipna
            const int64_t k = ins.first->second;
            const double y = v - comp[k];
            const double t = sum[k] + y;
            comp[k] = t - sum[k] - y;
            if (std::isnan(comp[k])) comp[k] = 0.0;  // +/-inf: keep the infinite sum
            sum[k] = t;
        });
    }
    pair_index = {};

    // Pairs sorted by (gene_a, gene_b) name, as groupby sorts its keys.
    const auto n_names = static_cast<int32_t>(names.size());
    std::vector<int32_t> order(n_names), rank(n_names);
    for (int32_t i = 0; i < n_names; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](int32_t x, int32_t y) { return names[x] < names[y]; });
    for (int32_t r = 0; r < n_names; ++r) rank[order[r]] = r;
    std::vector<int64_t> kept;
    for (int64_t k = 0; k < static_cast<int64_t>(pair_key.size()); ++k)
        if (std::fabs(sum[k]) >= cor_cutoff) kept.push_back(k);
    auto ranked = [&](int64_t k) {
        return (static_cast<uint64_t>(rank[pair_key[k] >> 32]) << 32) |
               static_cast<uint32_t>(rank[static_cast<uint32_t>(pair_key[k])]);
    };
    std::sort(kept.begin(), kept.end(),
              [&](int64_t x, int64_t y) { return ranked(x) < ranked(y); });

    MergedEdges res;
    std::vector<char> seen(n_names, 0);
    for (const int64_t k : kept) {
        const auto ia = static_cast<int32_t>(pair_key[k] >> 32);
        const auto ib = static_cast<int32_t>(static_cast<uint32_t>(pair_key[k]));
        const double s = sum[k];
        append_csv_field(res.text, names[ia]);
        res.text += '\t';
        append_csv_field(res.text, names[ib]);
        res.text += '\t';
        append_py_repr(res.text, s);
        res.text += s > 0 ? "\t1\n" : s < 0 ? "\t-1\n" : "\t0\n";
        res.node_count += !seen[ia];
        seen[ia] = 1;
        res.node_count += !seen[ib];
        seen[ib] = 1;
    }
    res.edge_count = static_cast<int64_t>(kept.size());
    return res;
}

namespace {

struct KeyedCor {
    uint64_t key;  // (first << 32) | second
    double cor;
};

// Interned gene names; the views in `ids` point into `names` (a deque, so
// they stay valid as names are added).
struct NameTable {
    std::deque<std::string> names;
    std::unordered_map<std::string_view, uint32_t> ids;

    uint32_t id_of(std::string_view s) {
        const auto it = ids.find(s);
        if (it != ids.end()) return it->second;
        const auto id = static_cast<uint32_t>(names.size());
        names.emplace_back(s);
        ids.emplace(names.back(), id);
        return id;
    }
};

// Number of lines in [b, e) (a last line without '\n' counts).
size_t count_lines(const char *b, const char *e) {
    const auto n = static_cast<size_t>(std::count(b, e, '\n'));
    return n + (e > b && e[-1] != '\n');
}

// Rows of the lines [b, e) as (local gene id, local gene id, correlation);
// line_no is the number of the line before b (for error messages).
void parse_diff_lines(const EdgeFileSpec &f, const char *b, const char *e, int64_t line_no,
                      NameTable &local, std::vector<KeyedCor> &rows) {
    const int want = f.cor_column + 1;
    std::vector<const char *> fs(want), fe(want);
    std::string sa, sb, sc;
    std::string last_a;  // edge lists are often grouped by the first gene
    uint32_t last_ia = 0;
    bool have_last = false;
    while (b < e) {
        const char *nl = static_cast<const char *>(std::memchr(b, '\n', e - b));
        const char *le = nl ? nl : e;
        const char *lb = b;
        b = nl ? nl + 1 : e;
        ++line_no;
        if (le > lb && le[-1] == '\r') --le;
        if (lb == le) continue;  // blank lines are skipped
        const int n = split_tabs(lb, le, fs.data(), fe.data(), want);
        const std::string_view a = unquote(fs[0], fe[0], sa);
        const std::string_view g = n > 1 ? unquote(fs[1], fe[1], sb) : std::string_view();
        if (n < 2 || is_na_token(a) || is_na_token(g))
            throw std::runtime_error(f.path + ": line " + std::to_string(line_no) +
                                     ": missing gene name");
        const double v = n > f.cor_column
                             ? parse_value(unquote(fs[f.cor_column], fe[f.cor_column], sc),
                                           f.path, line_no)
                             : std::numeric_limits<double>::quiet_NaN();
        if (!have_last || a != last_a) {
            last_ia = local.id_of(a);
            last_a.assign(a);
            have_last = true;
        }
        const uint32_t ib = local.id_of(g);
        rows.push_back({(static_cast<uint64_t>(last_ia) << 32) | ib, v});
    }
}

// The lines [b, e) parsed by up to n_threads threads, each on a run of whole
// lines with its own name table; the rows are appended to `rows` in file
// order with the ids of `names`. Line numbers start after line_no.
void parse_diff_block(const EdgeFileSpec &f, const char *b, const char *e, int64_t line_no,
                      NameTable &names, std::vector<KeyedCor> &rows, int n_threads,
                      size_t min_bytes) {
    using gcen_parallel::run_parallel;
    const size_t T = gcen_parallel::workers_for(static_cast<size_t>(e - b), n_threads,
                                                std::max<size_t>(1, min_bytes));
    std::vector<const char *> cut(T + 1, e);
    cut[0] = b;
    for (size_t t = 1; t < T; ++t) {
        const char *p = std::max(cut[t - 1], b + (e - b) * static_cast<ptrdiff_t>(t) /
                                                     static_cast<ptrdiff_t>(T));
        const char *nl = p < e ? static_cast<const char *>(std::memchr(p, '\n', e - p)) : nullptr;
        cut[t] = nl ? nl + 1 : e;
    }
    std::vector<int64_t> first_line(T + 1, line_no);
    std::vector<size_t> n_lines(T);
    run_parallel(T, [&](size_t t) { n_lines[t] = count_lines(cut[t], cut[t + 1]); });
    for (size_t t = 0; t < T; ++t) first_line[t + 1] = first_line[t] + n_lines[t];

    std::vector<NameTable> local(T);
    std::vector<std::vector<KeyedCor>> part(T);
    // the first error in file order is rethrown (run_parallel, by thread)
    run_parallel(T, [&](size_t t) {
        parse_diff_lines(f, cut[t], cut[t + 1], first_line[t], local[t], part[t]);
    });
    std::vector<std::vector<uint32_t>> remap(T);
    std::vector<size_t> offset(T + 1, rows.size());
    for (size_t t = 0; t < T; ++t) {  // in order, so ids do not depend on T
        remap[t].reserve(local[t].names.size());
        for (const auto &name : local[t].names) remap[t].push_back(names.id_of(name));
        offset[t + 1] = offset[t] + part[t].size();
    }
    rows.resize(offset[T]);
    run_parallel(T, [&](size_t t) {
        const auto &m = remap[t];
        KeyedCor *out = rows.data() + offset[t];
        for (const KeyedCor &r : part[t])
            *out++ = {(static_cast<uint64_t>(m[r.key >> 32]) << 32) |
                          m[static_cast<uint32_t>(r.key)],
                      r.cor};
        part[t] = {};
    });
}

// Rows of an edge-list file as (gene id, gene id, correlation), gene ids
// interned in `names`; the file is read in blocks of about kBlock bytes of
// whole lines, each parsed by up to n_threads threads (at least kBlock / 64
// bytes each).
void read_diff_rows(const EdgeFileSpec &f, NameTable &names, std::vector<KeyedCor> &rows,
                    int n_threads, size_t kBlock) {
    if (f.cor_column < 2) throw std::invalid_argument("cor_column must be >= 2");
    std::FILE *fh = std::fopen(f.path.c_str(), "rb");
    if (!fh) throw std::runtime_error("cannot open " + f.path);
    kBlock = std::max<size_t>(kBlock, 64);
    std::vector<char> buf;
    size_t have = 0;  // bytes of an unfinished line carried to the next block
    int64_t line_no = 0;
    bool header_pending = f.skip_header;
    try {
        for (;;) {
            buf.resize(have + kBlock);
            const size_t got = std::fread(buf.data() + have, 1, kBlock, fh);
            if (got == 0 && std::ferror(fh)) throw std::runtime_error("error reading " + f.path);
            const bool eof = got == 0;
            const char *b = buf.data();
            const char *end = b + have + got;
            const char *stop = end;  // parse [b, stop): whole lines
            if (!eof) {
                stop = b;
                for (const char *p = end; p > b; --p)
                    if (p[-1] == '\n') { stop = p; break; }
            }
            if (header_pending && b < stop) {
                const char *nl = static_cast<const char *>(std::memchr(b, '\n', stop - b));
                b = nl ? nl + 1 : stop;
                ++line_no;
                header_pending = false;
            }
            parse_diff_block(f, b, stop, line_no, names, rows, n_threads, kBlock / 64);
            line_no += static_cast<int64_t>(count_lines(b, stop));
            have = static_cast<size_t>(end - stop);
            std::memmove(buf.data(), stop, have);
            if (eof) break;
        }
    } catch (...) {
        std::fclose(fh);
        throw;
    }
    std::fclose(fh);
}

// Re-keys rows by sorted-name rank as (min << 32) | max, sorts them (stable)
// and keeps the last row of each edge.
void rank_and_dedup(std::vector<KeyedCor> &rows, const std::vector<uint32_t> &rank,
                    int n_threads) {
    const size_t n = rows.size();
    const size_t T = gcen_parallel::workers_for(n, n_threads);
    gcen_parallel::run_parallel(T, [&](size_t t) {
        for (size_t i = n * t / T; i < n * (t + 1) / T; ++i) {
            KeyedCor &r = rows[i];
            uint64_t x = rank[r.key >> 32], y = rank[static_cast<uint32_t>(r.key)];
            if (x > y) std::swap(x, y);
            r.key = (x << 32) | y;
        }
    });
    gcen_parallel::parallel_sort(
        rows, [](const KeyedCor &u, const KeyedCor &v) { return u.key < v.key; }, n_threads,
        /*stable=*/true);
    size_t out = 0;
    for (size_t i = 0; i < n; ++i) {
        if (i + 1 < n && rows[i + 1].key == rows[i].key) continue;
        rows[out++] = rows[i];
    }
    rows.resize(out);
    rows.shrink_to_fit();
}

}  // namespace

DiffTable diff_join_edge_files(const EdgeFileSpec &fa, const EdgeFileSpec &fb, int n_threads,
                               size_t block_bytes) {
    NameTable table;
    std::vector<KeyedCor> ra, rb;
    read_diff_rows(fa, table, ra, n_threads, block_bytes);
    read_diff_rows(fb, table, rb, n_threads, block_bytes);
    table.ids = {};
    std::vector<std::string> names(std::make_move_iterator(table.names.begin()),
                                   std::make_move_iterator(table.names.end()));
    table.names = {};

    const auto n_names = names.size();
    std::vector<uint32_t> order(n_names), rank(n_names);
    for (size_t i = 0; i < n_names; ++i) order[i] = static_cast<uint32_t>(i);
    std::sort(order.begin(), order.end(),
              [&](uint32_t x, uint32_t y) { return names[x] < names[y]; });
    for (size_t r = 0; r < n_names; ++r) rank[order[r]] = static_cast<uint32_t>(r);
    rank_and_dedup(ra, rank, n_threads);
    rank_and_dedup(rb, rank, n_threads);
    DiffTable t;
    t.names.reserve(n_names);
    for (const uint32_t i : order) t.names.push_back(std::move(names[i]));
    names = {};

    size_t cap = 0;  // rows of the join
    for (size_t i = 0, j = 0; i < ra.size() || j < rb.size(); ++cap) {
        if (j == rb.size() || (i < ra.size() && ra[i].key < rb[j].key)) ++i;
        else if (i == ra.size() || rb[j].key < ra[i].key) ++j;
        else { ++i; ++j; }
    }
    t.lo.reserve(cap); t.hi.reserve(cap);
    t.cor_a.reserve(cap); t.cor_b.reserve(cap); t.flags.reserve(cap);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    size_t i = 0, j = 0;
    while (i < ra.size() || j < rb.size()) {
        uint64_t key;
        double ca = nan, cb = nan;
        uint8_t fl = 0;
        if (j == rb.size() || (i < ra.size() && ra[i].key < rb[j].key)) {
            key = ra[i].key; ca = ra[i++].cor; fl = 1;
        } else if (i == ra.size() || rb[j].key < ra[i].key) {
            key = rb[j].key; cb = rb[j++].cor; fl = 2;
        } else {
            key = ra[i].key; ca = ra[i++].cor; cb = rb[j++].cor; fl = 3;
        }
        t.lo.push_back(static_cast<int32_t>(key >> 32));
        t.hi.push_back(static_cast<int32_t>(static_cast<uint32_t>(key)));
        t.cor_a.push_back(ca);
        t.cor_b.push_back(cb);
        t.flags.push_back(fl);
    }
    return t;
}

std::vector<int64_t> DiffTable::write(const std::string &path, const double *z, const double *p,
                                      const double *q, size_t n_shared, double pval,
                                      bool use_fdr, int n_threads) const {
    static const char *const kCategory[4] = {"A_only", "B_only", "shared_changed",
                                             "shared_stable"};
    constexpr size_t kBlockRows = size_t(1) << 15;
    const size_t n_rows = flags.size();
    const size_t n_blocks = (n_rows + kBlockRows - 1) / kBlockRows;
    // index of the first shared edge of each block into z / p / q
    std::vector<size_t> first_shared(n_blocks + 1, 0);
    for (size_t blk = 0; blk < n_blocks; ++blk) {
        size_t c = 0;
        for (size_t r = blk * kBlockRows; r < std::min(n_rows, (blk + 1) * kBlockRows); ++r)
            c += flags[r] == 3;
        first_shared[blk + 1] = first_shared[blk] + c;
    }
    if (first_shared[n_blocks] != n_shared)
        throw std::invalid_argument(std::to_string(n_shared) + " z / p values for " +
                                    std::to_string(first_shared[n_blocks]) + " shared edges");
    std::vector<std::string> qnames;
    qnames.reserve(names.size());
    for (const auto &s : names) {
        std::string f;
        append_csv_field(f, s);
        qnames.push_back(std::move(f));
    }
    std::vector<std::array<int64_t, 4>> block_counts(n_blocks, {0, 0, 0, 0});

    std::FILE *fh = std::fopen(path.c_str(), "wb");
    if (!fh) throw std::runtime_error("cannot write " + path);
    static const std::string header =
        "gene_a\tgene_b\tcor_a\tcor_b\tz_score\tp_value\tfdr\tcategory\n";
    bool ok = std::fwrite(header.data(), 1, header.size(), fh) == header.size();
    try {
        ok = ok && gcen_parallel::write_blocks(
                       fh, n_rows, n_threads,
                       [&](std::string &buf, size_t begin, size_t end) {
                           const size_t blk = begin / kBlockRows;
                           size_t k = first_shared[blk];
                           auto &counts = block_counts[blk];
                           for (size_t r = begin; r < end; ++r) {
                               buf += qnames[lo[r]];
                               buf += '\t';
                               buf += qnames[hi[r]];
                               buf += '\t';
                               append_py_repr(buf, cor_a[r]);
                               buf += '\t';
                               append_py_repr(buf, cor_b[r]);
                               buf += '\t';
                               int cat;
                               if (flags[r] == 3) {
                                   append_py_repr(buf, z[k]);
                                   buf += '\t';
                                   append_py_repr(buf, p[k]);
                                   buf += '\t';
                                   append_py_repr(buf, q[k]);
                                   cat = (use_fdr ? q[k] : p[k]) <= pval ? 2 : 3;  // NaN: stable
                                   ++k;
                               } else {
                                   buf += "\t\t";
                                   cat = flags[r] == 1 ? 0 : 1;
                               }
                               buf += '\t';
                               buf += kCategory[cat];
                               buf += '\n';
                               ++counts[cat];
                           }
                       },
                       kBlockRows);
    } catch (...) {
        std::fclose(fh);
        throw;
    }
    if (std::fclose(fh) != 0) ok = false;
    if (!ok) throw std::runtime_error("error writing " + path);
    std::vector<int64_t> counts(4, 0);
    for (const auto &c : block_counts)
        for (int i = 0; i < 4; ++i) counts[i] += c[i];
    return counts;
}

std::vector<double> bh_adjust(const double *p, size_t n, int n_threads) {
    struct PI {
        double p;
        size_t i;
    };
    std::vector<PI> v;
    v.reserve(n);
    for (size_t i = 0; i < n; ++i)
        if (std::isfinite(p[i])) v.push_back({p[i], i});
    gcen_parallel::parallel_sort(
        v, [](const PI &a, const PI &b) { return a.p < b.p || (a.p == b.p && a.i < b.i); },
        n_threads);
    std::vector<double> q(n, std::numeric_limits<double>::quiet_NaN());
    const auto m = static_cast<double>(v.size());
    double run = std::numeric_limits<double>::infinity();
    for (size_t k = v.size(); k-- > 0;) {
        run = std::min(run, v[k].p * m / static_cast<double>(k + 1));
        q[v[k].i] = std::min(run, 1.0);
    }
    return q;
}
