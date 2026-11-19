#include "table_writer.hpp"

#include <algorithm>
#include <cstdio>
#include <exception>
#include <stdexcept>

#include "edge_text.hpp"
#include "py_repr.hpp"
#include "parallel.hpp"

namespace {

// Validated name / label strings, CSV-quoted once instead of once per row.
std::vector<std::string> quoted(const std::vector<std::string> &v) {
    std::vector<std::string> out;
    out.reserve(v.size());
    for (const auto &s : v) {
        std::string q;
        edge_text::append_csv_field(q, s);
        out.push_back(std::move(q));
    }
    return out;
}

const std::string &at(const std::vector<std::string> &v, int64_t k) {
    if (k < 0 || static_cast<uint64_t>(k) >= v.size())
        throw std::out_of_range("code " + std::to_string(k) + " outside 0.." +
                                std::to_string(v.size()));
    return v[static_cast<size_t>(k)];
}

void format_rows(std::string &out, const std::vector<std::string> &names,
                 const std::vector<std::vector<std::string>> &labels, const int64_t *a,
                 const int64_t *b, const std::vector<TableColumn> &cols, size_t begin,
                 size_t end) {
    char ibuf[24];
    for (size_t r = begin; r < end; ++r) {
        out += at(names, a[r]);
        out += '\t';
        out += at(names, b[r]);
        for (size_t c = 0; c < cols.size(); ++c) {
            out += '\t';
            const TableColumn &col = cols[c];
            switch (col.kind) {
            case TableColumn::Float:
                append_py_repr(out, col.f[r]);
                break;
            case TableColumn::Int: {
                int n = std::snprintf(ibuf, sizeof ibuf, "%lld",
                                      static_cast<long long>(col.i[r]));
                out.append(ibuf, static_cast<size_t>(n));
                break;
            }
            case TableColumn::Label:
                out += at(labels[c], col.i[r]);
                break;
            }
        }
        out += '\n';
    }
}

}  // namespace

void write_coded_table(const std::string &path, const std::string &header,
                       const std::vector<std::string> &names, const int64_t *a,
                       const int64_t *b, size_t n_rows, const std::vector<TableColumn> &cols,
                       int n_threads) {
    const std::vector<std::string> qnames = quoted(names);
    std::vector<std::vector<std::string>> qlabels(cols.size());
    for (size_t c = 0; c < cols.size(); ++c)
        if (cols[c].kind == TableColumn::Label) qlabels[c] = quoted(*cols[c].labels);

    std::FILE *fh = std::fopen(path.c_str(), "wb");
    if (!fh) throw std::runtime_error("cannot write " + path);
    bool ok = std::fwrite(header.data(), 1, header.size(), fh) == header.size();

    try {
        ok = ok && gcen_parallel::write_blocks(
                       fh, n_rows, n_threads, [&](std::string &out, size_t begin, size_t end) {
                           format_rows(out, qnames, qlabels, a, b, cols, begin, end);
                       });
    } catch (...) {
        std::fclose(fh);
        throw;
    }
    if (std::fclose(fh) != 0) ok = false;
    if (!ok) throw std::runtime_error("error writing " + path);
}
