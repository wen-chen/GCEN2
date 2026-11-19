#include "layer_file.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string_view>
#include <unordered_map>
#if __has_include(<charconv>)
#include <charconv>
#endif

#include "edge_text.hpp"
#include "parse_double.hpp"

namespace {

struct Decline {
    std::string reason;
};

// [-]digits[.digits*] | [-].digits, then optional e[+-]digits
bool plain_decimal(std::string_view s) {
    size_t i = 0;
    if (i < s.size() && s[i] == '-') ++i;
    size_t int_digits = 0, frac_digits = 0;
    while (i < s.size() && s[i] >= '0' && s[i] <= '9') { ++i; ++int_digits; }
    if (i < s.size() && s[i] == '.') {
        ++i;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { ++i; ++frac_digits; }
    }
    if (int_digits + frac_digits == 0) return false;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        ++i;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) ++i;
        size_t exp_digits = 0;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { ++i; ++exp_digits; }
        if (exp_digits == 0) return false;
    }
    return i == s.size();
}

double parse_weight(std::string_view s, int64_t line_no) {
    if (!plain_decimal(s))
        throw Decline{"line " + std::to_string(line_no) + ": weight '" + std::string(s) +
                      "' is not a plain decimal number"};
    double v = 0.0;
    if (!parse_double(s.data(), s.data() + s.size(), v))
        throw Decline{"line " + std::to_string(line_no) + ": weight out of range"};
    return std::fabs(v);
}

}  // namespace

LayerFileResult read_layer_file(const LayerFileSpec &spec) {
    LayerFileResult r;
    if (spec.n_cols < 2 || spec.weight_col >= spec.n_cols) {
        r.reason = "bad column layout";
        return r;
    }
    std::deque<std::string> store;  // stable storage for the map keys
    std::unordered_map<std::string_view, int32_t> ids;
    ids.reserve(1 << 16);
    auto code = [&](std::string_view g) -> int32_t {
        auto it = ids.find(g);
        if (it != ids.end()) return it->second;
        if (store.size() >= static_cast<size_t>(INT32_MAX)) throw Decline{"too many genes"};
        store.emplace_back(g);
        const int32_t id = static_cast<int32_t>(store.size() - 1);
        ids.emplace(store.back(), id);
        return id;
    };

    const int max_fields = spec.n_cols + 1;  // one more, to see a longer row
    std::vector<const char *> starts(static_cast<size_t>(max_fields)),
        ends(static_cast<size_t>(max_fields));
    int64_t line_no = 0;
    bool header_pending = spec.skip_header;
    try {
        edge_text::for_each_line(spec.path, [&](const char *b, const char *e) {
            ++line_no;
            if (line_no == 1 && e - b >= 3 && std::memcmp(b, "\xEF\xBB\xBF", 3) == 0)
                throw Decline{"UTF-8 byte order mark"};
            if (b == e) return;  // pandas skips blank lines
            const size_t len = static_cast<size_t>(e - b);
            if (std::memchr(b, '"', len) || std::memchr(b, '\r', len))
                throw Decline{"line " + std::to_string(line_no) + ": quote or CR"};
            if (header_pending) {
                header_pending = false;
                return;
            }
            // fields split at sep (as split_tabs, for any one-char separator)
            int n = 0;
            const char *p = b;
            while (n < max_fields) {
                const char *t = static_cast<const char *>(std::memchr(p, spec.sep, e - p));
                starts[n] = p;
                ends[n] = t ? t : e;
                ++n;
                if (!t) break;
                p = t + 1;
            }
            if (n > spec.n_cols)
                throw Decline{"line " + std::to_string(line_no) + ": more fields than the first"};
            auto field = [&](int c) -> std::string_view {
                if (c >= n) return std::string_view();  // missing: NA, as pandas pads
                return std::string_view(starts[c], static_cast<size_t>(ends[c] - starts[c]));
            };
            const std::string_view ga = field(0), gb = field(1);
            if (edge_text::is_na_token(ga) || edge_text::is_na_token(gb))
                throw Decline{"line " + std::to_string(line_no) + ": missing gene name"};
            double w = 1.0;
            if (spec.weight_col >= 0) {
                const std::string_view ws = field(spec.weight_col);
                if (edge_text::is_na_token(ws))
                    throw Decline{"line " + std::to_string(line_no) + ": missing weight"};
                w = parse_weight(ws, line_no);
            }
            r.a.push_back(code(ga));
            r.b.push_back(code(gb));
            r.w.push_back(w);
        });
    } catch (const Decline &d) {
        r = LayerFileResult();
        r.reason = d.reason;
        return r;
    }
    r.names.assign(std::make_move_iterator(store.begin()), std::make_move_iterator(store.end()));
    r.ok = true;
    return r;
}
