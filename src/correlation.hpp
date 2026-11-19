#pragma once
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

/**
 * Compute all pairwise Pearson or Spearman correlations for an n_genes x n_samples matrix.
 *
 * Threading model: rows are distributed dynamically over the threads; output
 * is sorted by (gene_i, gene_j).
 * GIL must be released by the caller before entering this function.
 *
 * @param data      Row-major float64 matrix, shape [n_genes, n_samples].
 * @param n_genes   Number of genes (rows).
 * @param n_samples Number of samples (columns).
 * @param spearman  If true, precompute ranks and compute Spearman correlation.
 * @param n_threads Number of parallel threads.
 * @param min_abs_cor   Drop pairs with |r| < min_abs_cor (0 disables).
 * @param max_p         Drop pairs with p > max_p (+inf disables).
 * @param positive_only Drop pairs with r < 0.
 * @param row_begin, row_end  Only pairs (i, j) with row_begin <= i < row_end
 *                  and j > i are computed (row_end < 0: all rows). Lets the
 *                  caller stream the n(n-1)/2 pairs in blocks.
 * @return Vector of (gene_i, gene_j, correlation, p_value) tuples that pass
 *         the filters. Filtering happens inside the worker threads so that
 *         memory scales with the number of retained edges, not with
 *         n_genes^2 (needed for >10k genes, as in GCEN v1).
 */
std::vector<std::tuple<int, int, double, double>>
compute_pairwise_correlation(const double *data,
                             int n_genes,
                             int n_samples,
                             bool spearman,
                             int n_threads,
                             double min_abs_cor = 0.0,
                             double max_p = std::numeric_limits<double>::infinity(),
                             bool positive_only = false,
                             int row_begin = 0,
                             int row_end = -1);

/** One computed gene pair. */
struct CorEdge {
    int i, j;
    double corr, pval;
};

/**
 * The per-gene vectors the correlation kernel works on: ranks (Spearman,
 * as GCEN v1) or a copy of the values (Pearson). compute_pairwise_correlation
 * is prepare_rows + correlate_rows; a caller streaming many row blocks can
 * prepare once.
 */
std::vector<std::vector<double>> prepare_rows(const double *data, int n_genes, int n_samples,
                                              bool spearman);

/**
 * Threads kept alive across calls: run(job) executes job() on every thread
 * (the caller's included) and returns when all have finished. A caller
 * streaming many small row blocks pays the thread start-up once instead of
 * per block.
 */
class WorkerPool {
public:
    explicit WorkerPool(int n_threads);
    ~WorkerPool();
    WorkerPool(const WorkerPool &) = delete;
    WorkerPool &operator=(const WorkerPool &) = delete;
    int size() const { return static_cast<int>(threads_.size()) + 1; }
    void run(const std::function<void()> &job);

private:
    void loop();
    std::vector<std::thread> threads_;
    std::mutex mu_;
    std::condition_variable start_cv_, done_cv_;
    const std::function<void()> *job_ = nullptr;
    std::uint64_t generation_ = 0;
    int pending_ = 0;
    bool stop_ = false;
};

/** Pairs (i, j > i), row_begin <= i < row_end, passing the filters (see
 *  compute_pairwise_correlation), sorted by (i, j) unless `sorted` is false
 *  (then in no particular order). Uses `pool` when given (n_threads is then
 *  ignored), otherwise n_threads fresh threads. */
std::vector<CorEdge> correlate_rows(const std::vector<std::vector<double>> &matrix,
                                    int n_samples, int n_threads, double min_abs_cor,
                                    double max_p, bool positive_only, int row_begin,
                                    int row_end, WorkerPool *pool = nullptr,
                                    bool sorted = true);

/**
 * Fast path for many samples: every row centred and scaled to unit length
 * (z = (x - mean) / sqrt(sum (x - mean)^2), all zeros for a constant row), so
 * r(i, j) = z_i . z_j, computed in tiles by the TOM's product kernel
 * (gemm_nt, tom_gemm.hpp: fixed summation order, so the results do not
 * depend on the threads, the tiling or the instruction set). r agrees with
 * the per-pair kernel (GCEN v1's two-pass sum) to rounding (1e-15 to 1e-14), so
 * an edge whose |r| or p-value lies on a cutoff can differ.
 */
struct GemmRows {
    int n_genes = 0, n_samples = 0;
    std::vector<double> z;  // n_genes x n_samples, row-major
};

GemmRows prepare_gemm_rows(const std::vector<std::vector<double>> &rows, int n_samples);

/** As correlate_rows, on prepared GemmRows (r clipped to [-1, 1]). */
std::vector<CorEdge> correlate_rows_gemm(const GemmRows &rows, int n_threads,
                                         double min_abs_cor, double max_p, bool positive_only,
                                         int row_begin, int row_end, WorkerPool *pool = nullptr,
                                         bool sorted = true);

/** Samples from which kernel "auto" uses the product kernel. */
constexpr int kGemmMinSamples = 64;

/** Whether `kernel` ("auto", "pairwise", "gemm") uses the product kernel for
 *  n_samples samples; throws std::invalid_argument for another name. */
bool use_gemm_kernel(const std::string &kernel, int n_samples);

/**
 * Dropout-aware Pearson: per gene pair, the correlation over a subset of the
 * samples, with the t-test p-value on the size of that subset. Pairs with
 * fewer than 3 such samples, or constant on them, get r = 0 and p = 1.
 *
 *  - BothNonzero ("pearson_nonzero"): the samples where both genes are
 *    non-zero (zeros treated as missing values).
 *  - AnyNonzero ("pearson_pairwise_complete"): the samples where at least one
 *    gene is non-zero (only double zeros dropped). Biased towards negative
 *    correlations on sparse data: the samples left with one gene at 0 and
 *    the other non-zero are anti-correlated. Pairs without double zeros give
 *    exactly the plain Pearson result.
 */
enum class NonzeroMode { BothNonzero, AnyNonzero };

struct NonzeroRows {
    NonzeroMode mode = NonzeroMode::BothNonzero;
    int n_samples = 0;
    // AnyNonzero: x - mean(x) over all samples; BothNonzero: x - mean of the
    // non-zero values, NaN where x is 0
    std::vector<std::vector<double>> values;
    std::vector<double> mean;
    std::vector<double> ss;                     // AnyNonzero: sum (x - mean)^2
    std::vector<std::vector<uint64_t>> zero;    // AnyNonzero: bit k: sample k is 0
    std::vector<std::vector<int>> nz;           // BothNonzero: non-zero samples
};

NonzeroRows prepare_nonzero_rows(const double *data, int n_genes, int n_samples,
                                 NonzeroMode mode);

// r and the number of samples used for genes i, j.
void nonzero_pearson(const NonzeroRows &rows, int i, int j, double &corr, double &n_eff);

std::vector<CorEdge> correlate_rows_nonzero(const NonzeroRows &rows, int n_threads,
                                            double min_abs_cor, double max_p,
                                            bool positive_only, int row_begin, int row_end,
                                            WorkerPool *pool = nullptr,
                                            bool sorted = true);
