#pragma once
/**
 * Tab-separated table writer for edge tables whose genes are integer codes
 * into a name list, written byte for byte as pandas
 * ``DataFrame.to_csv(path, sep="\t", index=False)`` writes the decoded table:
 * strings quoted as QUOTE_MINIMAL (append_csv_field), floats as Python
 * repr() (append_py_repr; NaN as an empty field), integers in decimal.
 *
 * Rows are formatted by several threads in blocks and written in order.
 * Used by gcen2.multilayer.fuse_to_file.
 */
#include <cstdint>
#include <string>
#include <vector>

struct TableColumn {
    enum Kind { Float, Int, Label };
    Kind kind = Float;
    const double *f = nullptr;    // Float: one value per row
    const int64_t *i = nullptr;   // Int: value; Label: index into labels
    const std::vector<std::string> *labels = nullptr;  // Label only
};

// ``header`` is written as given (with its newline). Gene columns come first:
// names[a[r]], names[b[r]]. Throws std::runtime_error when the file cannot be
// written and std::out_of_range for a code outside names / labels.
void write_coded_table(const std::string &path, const std::string &header,
                       const std::vector<std::string> &names, const int64_t *a,
                       const int64_t *b, size_t n_rows, const std::vector<TableColumn> &cols,
                       int n_threads);
