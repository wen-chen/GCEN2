#pragma once
/**
 * Edge-list file to gene codes, for gcen2.multilayer: the genes of the two
 * gene columns as int32 codes into a list of names (first-seen order) and
 * the |weight| per row, read as gcen2.multilayer.read_layer (pandas
 * ``read_csv(sep=sep, header=None, dtype=str)`` and ``pd.to_numeric``)
 * reads them, without building strings per row.
 *
 * The header line and the column layout are decided by the caller (from
 * the first record, read by pandas). Only inputs whose pandas reading is
 * certain are handled; otherwise ok = false and the caller falls back to
 * pandas: a UTF-8 BOM, a '"' or a '\r' inside a line, a row longer than
 * the first, a missing (NA) gene or weight, a weight that is not a plain
 * decimal number ([-]digits[.digits][e[+-]digits]; pandas parses a leading
 * '+' differently) or that overflows.
 */
#include <cstdint>
#include <string>
#include <vector>

struct LayerFileSpec {
    std::string path;
    char sep = '\t';
    bool skip_header = false;  // the first non-blank line is a header
    int n_cols = 2;            // fields of the first record
    int weight_col = -1;       // -1: weight 1 for every row
};

struct LayerFileResult {
    bool ok = false;
    std::string reason;  // why the file was not handled (ok == false)
    std::vector<std::string> names;
    std::vector<int32_t> a, b;
    std::vector<double> w;
};

LayerFileResult read_layer_file(const LayerFileSpec &spec);
