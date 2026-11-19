#pragma once
// Text helpers shared by the edge-list file readers and writers
// (edge_files.cpp, shuffle_file.cpp): lines, tab-separated fields, pandas'
// unquoting and missing-value tokens, to_csv field quoting.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace edge_text {

// Calls fn(begin, end) for every line of the file ('\n' and a trailing '\r'
// removed), reading in large chunks.
template <typename Fn>
inline void for_each_line(const std::string &path, Fn &&fn) {
    std::FILE *fh = std::fopen(path.c_str(), "rb");
    if (!fh) throw std::runtime_error("cannot open " + path);
    constexpr size_t kChunk = 1 << 22;
    std::vector<char> buf(kChunk);
    std::string carry;
    auto emit = [&](const char *b, const char *e) {
        if (e > b && *(e - 1) == '\r') --e;
        fn(b, e);
    };
    try {
        for (;;) {
            const size_t got = std::fread(buf.data(), 1, kChunk, fh);
            if (got == 0) break;
            const char *p = buf.data();
            const char *end = p + got;
            while (p < end) {
                const char *nl = static_cast<const char *>(std::memchr(p, '\n', end - p));
                if (!nl) { carry.append(p, end); break; }
                if (!carry.empty()) {
                    carry.append(p, nl);
                    emit(carry.data(), carry.data() + carry.size());
                    carry.clear();
                } else {
                    emit(p, nl);
                }
                p = nl + 1;
            }
        }
        if (!carry.empty()) emit(carry.data(), carry.data() + carry.size());
    } catch (...) {
        std::fclose(fh);
        throw;
    }
    const bool failed = std::ferror(fh) != 0;
    std::fclose(fh);
    if (failed) throw std::runtime_error("error reading " + path);
}

// A field as pandas' C tokenizer reads it: a field starting with '"' is
// quoted; inside the quotes "" is a '"', and after the closing quote the
// rest of the field is taken literally ("q"t" -> qt"). A field not starting
// with '"' is taken as is. (Lines are split at tabs and newlines first, so
// a quoted field holding a tab or newline is not supported.) Returns a view
// of the field itself when it is not quoted (no copy).
inline std::string_view unquote(const char *b, const char *e, std::string &scratch) {
    if (b == e || *b != '"') return std::string_view(b, static_cast<size_t>(e - b));
    scratch.clear();
    bool quoted = true;
    for (const char *p = b + 1; p < e; ++p) {
        if (quoted && *p == '"') {
            if (p + 1 < e && *(p + 1) == '"') {
                scratch.push_back('"');
                ++p;
            } else {
                quoted = false;  // closing quote
            }
        } else {
            scratch.push_back(*p);
        }
    }
    return scratch;
}

// Splits [b, e) at tabs into at most `max_fields` field ranges.
inline int split_tabs(const char *b, const char *e, const char **starts, const char **ends,
               int max_fields) {
    int n = 0;
    const char *p = b;
    while (n < max_fields) {
        const char *t = static_cast<const char *>(std::memchr(p, '\t', e - p));
        starts[n] = p;
        ends[n] = t ? t : e;
        ++n;
        if (!t) break;
        p = t + 1;
    }
    return n;
}

// pandas' default missing-value strings.
inline bool is_na_token(std::string_view s) {
    static const std::unordered_set<std::string_view> na = {
        "",     "#N/A", "#N/A N/A", "#NA",  "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND",
        "1.#QNAN", "<NA>", "N/A",  "NA",   "NULL",    "NaN",      "None", "n/a",  "nan",
        "null"};
    return na.count(s) != 0;
}

inline double parse_value(std::string_view s, const std::string &path, int64_t line_no) {
    if (is_na_token(s)) return std::numeric_limits<double>::quiet_NaN();
    std::string tmp(s);
    char *end = nullptr;
    const double v = std::strtod(tmp.c_str(), &end);
    if (end == tmp.c_str() || *end != '\0')
        throw std::runtime_error(path + ": line " + std::to_string(line_no) +
                                 ": correlation '" + tmp + "' is not a number");
    return v;
}

// A field as pandas' to_csv writes it (csv.QUOTE_MINIMAL): quoted, with
// '"' doubled, when it holds a quote, tab or line break.
inline void append_csv_field(std::string &out, std::string_view f) {
    if (f.find_first_of("\"\t\r\n") == std::string::npos) {
        out += f;
        return;
    }
    out += '"';
    for (const char c : f) {
        if (c == '"') out += '"';
        out += c;
    }
    out += '"';
}

}  // namespace edge_text
