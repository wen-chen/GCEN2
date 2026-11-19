#pragma once
// Python repr() of a double, shared by the formatters that must match pandas
// output (pandas writes float columns with repr()).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#if __has_include(<charconv>)
#include <charconv>
#endif

// Python repr() of a double (shortest round-trip digits; fixed notation for
// decimal exponents -4 <= e < 16, else d.ddde+XX), as pandas writes floats.
inline void append_py_repr(std::string &s, double v) {
    if (std::isnan(v)) return;  // pandas writes NaN as an empty field
    if (std::isinf(v)) { s += v > 0 ? "inf" : "-inf"; return; }
    char buf[64];
    int len = 0;
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
    auto res = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::scientific);
    len = static_cast<int>(res.ptr - buf);
#else
    for (int prec = 0; prec <= 17; ++prec) {
        len = std::snprintf(buf, sizeof buf, "%.*e", prec, v);
        if (std::strtod(buf, nullptr) == v) break;
    }
#endif
    // buf = [-]d[.ddd]e(+|-)XX
    std::string m(buf, static_cast<size_t>(len));
    size_t epos = m.find('e');
    int exp10 = std::atoi(m.c_str() + epos + 1);
    std::string digits;
    bool neg = false;
    for (size_t i = 0; i < epos; ++i) {
        if (m[i] == '-') neg = true;
        else if (m[i] != '.') digits += m[i];
    }
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    if (neg) s += '-';
    int nd = static_cast<int>(digits.size());
    if (exp10 >= -4 && exp10 < 16) {
        if (exp10 < 0) {
            s += "0.";
            s.append(static_cast<size_t>(-exp10 - 1), '0');
            s += digits;
        } else if (nd <= exp10 + 1) {
            s += digits;
            s.append(static_cast<size_t>(exp10 + 1 - nd), '0');
            s += ".0";
        } else {
            s.append(digits, 0, static_cast<size_t>(exp10 + 1));
            s += '.';
            s.append(digits, static_cast<size_t>(exp10 + 1), std::string::npos);
        }
    } else {
        s += digits[0];
        if (nd > 1) { s += '.'; s.append(digits, 1, std::string::npos); }
        char eb[16];
        std::snprintf(eb, sizeof eb, "e%c%02d", exp10 < 0 ? '-' : '+', std::abs(exp10));
        s += eb;
    }
}

