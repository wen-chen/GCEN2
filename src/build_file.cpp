#include "build_file.hpp"

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "correlation.hpp"
#include "parse_double.hpp"
#include "network_stream.hpp"

namespace {

// Python's str.isspace() for ASCII characters (str.strip() removes these).
inline bool py_ascii_space(unsigned char c) {
    return c == ' ' || (c >= 9 && c <= 13) || (c >= 28 && c <= 31);
}

bool valid_utf8(const unsigned char *p, const unsigned char *e) {
    while (p < e) {
        const unsigned char c = *p;
        if (c < 0x80) { ++p; continue; }
        int len;
        uint32_t cp;
        if ((c & 0xE0) == 0xC0) { len = 2; cp = c & 0x1F; }
        else if ((c & 0xF0) == 0xE0) { len = 3; cp = c & 0x0F; }
        else if ((c & 0xF8) == 0xF0) { len = 4; cp = c & 0x07; }
        else return false;
        if (e - p < len) return false;
        for (int k = 1; k < len; ++k) {
            if ((p[k] & 0xC0) != 0x80) return false;
            cp = (cp << 6) | (p[k] & 0x3F);
        }
        // overlong forms, surrogates and values above U+10FFFF are invalid
        if ((len == 2 && cp < 0x80) || (len == 3 && cp < 0x800) || (len == 4 && cp < 0x10000) ||
            (cp >= 0xD800 && cp <= 0xDFFF) || cp > 0x10FFFF)
            return false;
        p += len;
    }
    return true;
}

// A plain decimal number, parsed exactly (correctly rounded, like Python and
// numpy). Anything else -> false.
bool parse_plain_number(const char *b, const char *e, double &v) {
    if (b == e) return false;
    const char *p = b;
    if (*p == '-') ++p;
    if (p == e || !((*p >= '0' && *p <= '9') || *p == '.')) return false;
    for (const char *q = p; q < e; ++q) {
        const char c = *q;
        if (!((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-'))
            return false;
    }
    return parse_double(b, e, v) && std::isfinite(v);
}

struct PlainMatrix {
    std::vector<std::string> genes;
    std::vector<double> values;  // genes x width, row-major
    int width = 0;
};

// read_matrix (gcen2._matrix_io) on the plain subset; false if outside it.
bool read_plain_matrix(const std::string &path, bool has_header, PlainMatrix &out) {
    std::FILE *fh = std::fopen(path.c_str(), "rb");
    if (!fh) return false;
    std::vector<double> values;
    std::vector<int32_t> lens;
    int64_t header_fields = -1;
    bool header_pending = has_header;
    bool ok = true;

    auto handle_line = [&](const char *b, const char *e) {
        if (!valid_utf8(reinterpret_cast<const unsigned char *>(b),
                        reinterpret_cast<const unsigned char *>(e))) {
            ok = false;
            return;
        }
        while (e > b && (*(e - 1) == '\r' || *(e - 1) == '\n')) --e;  // rstrip("\r\n")
        if (header_pending) {  // the first line, even if blank
            header_pending = false;
            header_fields = 1 + std::count(b, e, '\t');
            return;
        }
        bool blank = true, ascii_text = false;
        for (const char *p = b; p < e; ++p) {
            const auto c = static_cast<unsigned char>(*p);
            if (!py_ascii_space(c)) {
                blank = false;
                if (c < 0x80) { ascii_text = true; break; }
            }
        }
        if (blank) return;
        // non-ASCII whitespace could make Python treat the line as blank
        if (!ascii_text) { ok = false; return; }
        const char *t = static_cast<const char *>(std::memchr(b, '\t', e - b));
        const char *name_end = t ? t : e;
        if (name_end - b >= 2 && *b == '"' && *(name_end - 1) == '"')
            out.genes.emplace_back(b + 1, name_end - 1);
        else
            out.genes.emplace_back(b, name_end);
        int32_t n = 0;
        while (t) {
            const char *f = t + 1;
            t = static_cast<const char *>(std::memchr(f, '\t', e - f));
            double v;
            if (!parse_plain_number(f, t ? t : e, v)) { ok = false; return; }
            values.push_back(v);
            ++n;
        }
        lens.push_back(n);
    };

    constexpr size_t kChunk = 1 << 22;
    std::vector<char> buf(kChunk);
    std::string carry;
    bool any = false;
    while (ok) {
        const size_t got = std::fread(buf.data(), 1, kChunk, fh);
        if (got == 0) break;
        any = true;
        const char *p = buf.data();
        const char *end = p + got;
        while (ok && p < end) {
            const char *nl = static_cast<const char *>(std::memchr(p, '\n', end - p));
            if (!nl) { carry.append(p, end); break; }
            if (!carry.empty()) {
                carry.append(p, nl);
                handle_line(carry.data(), carry.data() + carry.size());
                carry.clear();
            } else {
                handle_line(p, nl);
            }
            p = nl + 1;
        }
    }
    const bool failed = std::ferror(fh) != 0;
    std::fclose(fh);
    if (failed || !ok) return false;
    if (!carry.empty()) {
        handle_line(carry.data(), carry.data() + carry.size());
        if (!ok) return false;
    }
    if (header_pending && !any) header_fields = 1;  // next(lines, "") -> one empty field

    int64_t width = 0;
    for (const int32_t n : lens) width = std::max<int64_t>(width, n);
    if (header_fields >= 0) {
        width = std::max<int64_t>(width, header_fields - 1);
        if (header_fields != width && header_fields - 1 != width) return false;  // Python reports
    }
    // build_to_file needs >= 3 samples; the kernels take int sizes
    if (out.genes.empty() || width < 3 || width > (1 << 30) ||
        out.genes.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        return false;
    out.width = static_cast<int>(width);
    bool ragged = false;
    for (const int32_t n : lens) ragged |= (n != width);
    if (!ragged) {
        out.values = std::move(values);
    } else {  // pad short rows with NaN (np.full(..., np.nan))
        out.values.assign(out.genes.size() * static_cast<size_t>(width),
                          std::numeric_limits<double>::quiet_NaN());
        size_t src = 0;
        for (size_t i = 0; i < lens.size(); ++i) {
            std::copy(values.begin() + src, values.begin() + src + lens[i],
                      out.values.begin() + i * width);
            src += lens[i];
        }
    }
    return true;
}

}  // namespace

BuildFileResult build_network_file(const std::string &in_path, const std::string &out_path,
                                   bool has_header, bool spearman, double cor, double pval,
                                   bool positive_only, bool all_edges, int n_threads,
                                   int64_t max_pairs_per_block, bool fdr, const TomSpec &tom,
                                   const std::string &kernel) {
    BuildFileResult res;
    PlainMatrix m;
    if (!read_plain_matrix(in_path, has_header, m)) return res;
    const int n = static_cast<int>(m.genes.size());
    WorkerPool pool(n_threads);  // started once, reused by every block
    const std::string method = spearman ? "spearman" : "pearson";
    const Correlate correlate =
        make_correlate(method, m.values.data(), n, m.width, pool, kernel);
    std::unique_ptr<TomScorer> scorer =
        make_tom_scorer(tom, method, m.values.data(), n, m.width, positive_only, pool);
    std::vector<double>().swap(m.values);

    std::FILE *fh = std::fopen(out_path.c_str(), "wb");
    if (!fh) return res;
    res.supported = true;
    EdgeFileWriter writer(fh, out_path, m.genes, fdr, pool, scorer != nullptr);
    NetworkOptions opt;
    opt.cor = cor;
    opt.pval = pval;
    opt.positive_only = positive_only;
    opt.fdr = fdr;
    opt.all_edges = all_edges;
    opt.max_pairs_per_block = max_pairs_per_block;
    opt.tom = scorer.get();
    opt.tom_min = tom.tom_min;
    stream_network(n, opt, pool.size(), correlate,
                   [&](const std::vector<CorEdge> &e, const std::vector<double> *q,
                       const std::vector<double> *t) { writer.write(e, q, t); });
    writer.close();
    if (scorer) res.tom = tom_info(*scorer);
    res.node_count = writer.node_count();
    res.edge_count = writer.edge_count();
    return res;
}
