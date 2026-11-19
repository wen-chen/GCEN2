#pragma once
// Correctly rounded text -> double, for strings already checked to be plain
// decimal numbers. std::from_chars where the standard library implements it
// for double (libstdc++ >= 11, MSVC); strtod elsewhere: Apple's libc++ declares
// the floating-point from_chars deleted. Both round correctly, so the result
// does not depend on the platform. (strtod reads the C locale's decimal point,
// which a Python process keeps as '.'.)

#include <cerrno>
#include <cmath>
#include <charconv>
#include <cstdlib>
#include <string>
#include <system_error>

inline bool parse_double(const char *b, const char *e, double &v) {
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L && !defined(_LIBCPP_VERSION)
    const auto res = std::from_chars(b, e, v);
    return res.ec == std::errc() && res.ptr == e;
#else
    const std::string tmp(b, e);
    char *end = nullptr;
    errno = 0;
    v = std::strtod(tmp.c_str(), &end);
    if (end != tmp.c_str() + tmp.size()) return false;
    // ERANGE also flags subnormal results, which from_chars accepts; fail only
    // on overflow and on underflow to zero, as from_chars does
    return !(errno == ERANGE && (std::isinf(v) || v == 0.0));
#endif
}
