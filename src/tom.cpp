#include "tom.hpp"

// The products go through gemm_nt (tom_gemm.hpp: fixed summation order, so
// values do not depend on the instruction set); tiles are spread over the
// WorkerPool, each product single-threaded, so values do not depend on the
// number of threads either.
#include "tom_gemm.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace {


constexpr int kCorRows = 64;    // rows per correlation tile
constexpr int kTileRows = 128;  // l_ij tiles: rows x columns
constexpr int kTileCols = 256;

// job(t) for t in [0, n) on the pool's threads.
template <class Job>
void parallel_for(WorkerPool &pool, int64_t n, Job job) {
    std::atomic<int64_t> next{0};
    pool.run([&] {
        for (int64_t t; (t = next.fetch_add(1)) < n;) job(t);
    });
}

// R's mean(): long double sum / n, then one correction pass.
double r_mean(const double *x, size_t n) {
    long double s = 0.0L;
    for (size_t i = 0; i < n; ++i) s += x[i];
    s /= static_cast<long double>(n);
    long double t = 0.0L;
    for (size_t i = 0; i < n; ++i) t += x[i] - s;
    s += t / static_cast<long double>(n);
    return static_cast<double>(s);
}

// R's seq.int(from, to, length.out = n) (C do_seq: symmetric from both ends).
std::vector<double> r_seq_int(double from, double to, int n) {
    std::vector<double> v(n);
    v[0] = from;
    if (n > 1) v[n - 1] = to;
    if (n > 2) {
        const double by = (to - from) / static_cast<double>(n - 1);
        for (int i = 1; i < n - 1; ++i)
            v[i] = i < n / 2 ? from + static_cast<double>(i) * by
                             : to - static_cast<double>(n - 1 - i) * by;
    }
    return v;
}

}  // namespace

std::vector<double> default_soft_powers() {
    std::vector<double> p;
    for (int i = 1; i <= 10; ++i) p.push_back(i);
    for (int i = 12; i <= 20; i += 2) p.push_back(i);
    return p;
}

double faq_power(int n_samples, bool signed_net) {
    const double u = n_samples < 20 ? 9 : n_samples < 30 ? 8 : n_samples < 40 ? 7 : 6;
    return signed_net ? 2 * u : u;
}

SoftThresholdFit scale_free_fit(const std::vector<double> &k, double power, int n_breaks) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SoftThresholdFit f{power, nan, nan, nan, nan, nan};
    const size_t n = k.size();
    if (n == 0) return f;
    f.mean_k = r_mean(k.data(), n);
    std::vector<double> s(k);
    std::sort(s.begin(), s.end());
    f.max_k = s.back();
    f.median_k = n % 2 ? s[n / 2] : (s[n / 2 - 1] + s[n / 2]) / 2.0;
    const double mn = s.front(), mx = s.back();
    if (!(mx > mn)) return f;  // hist() would stop: no fit

    // cut(k, n_breaks): equal-width breaks, outer ones widened by range/1000,
    // intervals (a, b]
    const int nb = n_breaks + 1;
    const double dx = mx - mn;
    std::vector<double> br = r_seq_int(mn, mx, nb);
    br[0] = mn - dx / 1000;
    br[nb - 1] = mx + dx / 1000;
    std::vector<std::vector<double>> bins(n_breaks);
    for (const double x : k) {
        int lo = 0, hi = nb - 1;  // .bincode, right = TRUE
        while (hi - lo >= 2) {
            const int mid = (hi + lo) / 2;
            if (x > br[mid]) lo = mid;
            else hi = mid;
        }
        bins[lo].push_back(x);
    }
    // hist(k, seq(min, max, length = n_breaks + 1))$mids for empty bins
    std::vector<double> b1(nb);
    b1[0] = mn;
    b1[nb - 1] = mx;
    for (int i = 1; i < nb - 1; ++i) b1[i] = mn + i * (dx / n_breaks);
    std::vector<double> xs(n_breaks), ys(n_breaks);
    for (int b = 0; b < n_breaks; ++b) {
        const double mid = 0.5 * (b1[b + 1] + b1[b]);
        double dk = bins[b].empty() ? mid : r_mean(bins[b].data(), bins[b].size());
        if (dk == 0) dk = mid;
        const double p = static_cast<double>(bins[b].size()) / static_cast<double>(n);
        xs[b] = std::log10(dk);
        ys[b] = std::log10(p + 1e-09);
    }
    // lm(log.p.dk ~ log.dk): slope and R^2
    double mxs = 0, mys = 0;
    for (int b = 0; b < n_breaks; ++b) { mxs += xs[b]; mys += ys[b]; }
    mxs /= n_breaks;
    mys /= n_breaks;
    double sxx = 0, sxy = 0, syy = 0;
    for (int b = 0; b < n_breaks; ++b) {
        const double u = xs[b] - mxs, v = ys[b] - mys;
        sxx += u * u;
        sxy += u * v;
        syy += v * v;
    }
    f.slope = sxy / sxx;
    f.r2 = sxy * sxy / (sxx * syy);
    return f;
}

// The adjacency in precision F (float for edge scores, double for modules)
// and its l tiles.
template <class F>
struct TomCore {
    int n = 0;
    std::vector<F> own;          // storage of A unless the caller gives it
    F *A = nullptr;              // r, then the adjacency (n x n, diagonal 0)
    std::vector<double> k;       // connectivity
    WorkerPool &pool;
    // l tiles of the current row tile, by column tile
    int cur_row_tile = -1;
    std::map<int, std::vector<F>> tiles;

    explicit TomCore(WorkerPool &p) : pool(p) {}

    // A = r (both triangles, NaN for genes without a valid row), in `into`
    // (n x n) when given
    // center = false: rows already centred (bicor rows), not re-centred
    void correlations(const std::vector<std::vector<double>> &rows, int n_samples,
                      F *into = nullptr, bool center = true) {
        n = static_cast<int>(rows.size());
        const int s = n_samples;
        std::vector<double> Z(static_cast<size_t>(n) * s);  // centred, unit-norm rows
        std::vector<uint8_t> ok(n, 1);
        for (int i = 0; i < n; ++i) {
            const double *x = rows[i].data();
            double mean = 0;
            if (center) {
                for (int t = 0; t < s; ++t) mean += x[t];
                mean /= s;
            }
            double ss = 0;
            for (int t = 0; t < s; ++t) ss += (x[t] - mean) * (x[t] - mean);
            if (!(ss > 0) || !std::isfinite(ss)) {
                ok[i] = 0;
                std::fill(&Z[static_cast<size_t>(i) * s], &Z[static_cast<size_t>(i + 1) * s], 0.0);
                continue;
            }
            const double inv = 1.0 / std::sqrt(ss);
            for (int t = 0; t < s; ++t) Z[static_cast<size_t>(i) * s + t] = (x[t] - mean) * inv;
        }
        if (into) {
            A = into;
        } else {
            own.assign(static_cast<size_t>(n) * n, F(0));
            A = own.data();
        }
        for (int64_t i = 0; i < n; ++i) A[i * n + i] = F(0);
        const F fnan = std::numeric_limits<F>::quiet_NaN();
        const int n_tiles = (n + kCorRows - 1) / kCorRows;
        parallel_for(pool, n_tiles, [&](int64_t t) {
            const int r0 = static_cast<int>(t) * kCorRows;
            const int R = std::min(kCorRows, n - r0);
            const int W = n - r0;
            std::vector<double> C(static_cast<size_t>(R) * W);
            gemm_nt(Z.data(), s, r0, R, r0, W, s, C.data());
            for (int a = 0; a < R; ++a) {
                const int i = r0 + a;
                for (int j = i + 1; j < n; ++j) {
                    F v = fnan;
                    if (ok[i] && ok[j])
                        v = static_cast<F>(std::max(
                            -1.0, std::min(1.0, C[static_cast<size_t>(a) * W + (j - r0)])));
                    A[static_cast<size_t>(i) * n + j] = v;
                    A[static_cast<size_t>(j) * n + i] = v;
                }
            }
        });
    }

    static double base(F r, bool signed_net) {
        const double x = r;
        return signed_net ? (1.0 + x) / 2.0 : std::fabs(x);
    }

    // connectivity for each power (sorted ascending), from A = r
    std::vector<std::vector<double>> connectivity(const std::vector<double> &powers,
                                                  bool signed_net) {
        const size_t P = powers.size();
        std::vector<double> steps(P);
        for (size_t p = 0; p < P; ++p) steps[p] = powers[p] - (p ? powers[p - 1] : 0.0);
        std::vector<std::vector<double>> kp(P, std::vector<double>(n, 0.0));
        parallel_for(pool, n, [&](int64_t i) {
            std::vector<double> acc(P, 0.0);
            const F *row = &A[static_cast<size_t>(i) * n];
            for (int j = 0; j < n; ++j) {
                if (j == i || std::isnan(row[j])) continue;
                const double b = base(row[j], signed_net);
                double cur = 1.0;
                for (size_t p = 0; p < P; ++p) {  // corxPrev * corx^step
                    const double st = steps[p];
                    cur *= st == 1.0 ? b : st == 2.0 ? b * b : std::pow(b, st);
                    acc[p] += cur;
                }
            }
            for (size_t p = 0; p < P; ++p) kp[p][i] = acc[p];
        });
        return kp;
    }

    void adjacency(double power, bool signed_net) {
        k.assign(n, 0.0);
        parallel_for(pool, n, [&](int64_t i) {
            F *row = &A[static_cast<size_t>(i) * n];
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                F a = 0;
                if (j != i && !std::isnan(row[j]))
                    a = static_cast<F>(std::pow(base(row[j], signed_net), power));
                row[j] = a;
                sum += a;
            }
            k[i] = sum;
        });
    }

    void ensure_tiles(int rt, const std::vector<int> &cts) {
        if (rt != cur_row_tile) {
            tiles.clear();
            cur_row_tile = rt;
        }
        std::vector<int> todo;
        for (const int c : cts)
            if (!tiles.count(c)) todo.push_back(c);
        if (todo.empty()) return;
        std::vector<std::vector<F>> out(todo.size());
        const int r0 = rt * kTileRows;
        parallel_for(pool, static_cast<int64_t>(todo.size()), [&](int64_t t) {
            out[t] = tile(r0, todo[t] * kTileCols);
        });
        for (size_t t = 0; t < todo.size(); ++t) tiles[todo[t]] = std::move(out[t]);
    }

    // l of rows [r0, r0 + 128) x columns [c0, c0 + 256) (one product, single-threaded)
    std::vector<F> tile(int r0, int c0) const {
        const int R = std::min(kTileRows, n - r0), C = std::min(kTileCols, n - c0);
        std::vector<F> out(static_cast<size_t>(R) * C);
        gemm_nt(A, n, r0, R, c0, C, n, out.data());
        return out;
    }

    // power > 0 as given, else picked from `powers` (pickSoftThreshold)
    TomInfo choose_power(double power, const std::vector<double> &powers, int n_samples,
                         bool signed_net) {
        TomInfo info;
        if (power > 0) {
            info.power = power;
            return info;
        }
        std::vector<double> pw(powers);
        std::sort(pw.begin(), pw.end());
        if (pw.empty() || pw.front() <= 0) throw std::invalid_argument("powers must be positive.");
        const auto kp = connectivity(pw, signed_net);
        for (size_t p = 0; p < pw.size(); ++p) {
            info.fits.push_back(scale_free_fit(kp[p], pw[p]));
            if (!info.estimated && info.fits.back().r2 > 0.85) {
                info.power = pw[p];
                info.estimated = true;
            }
        }
        if (!info.estimated) info.power = faq_power(n_samples, signed_net);
        return info;
    }

    // 1 - TOM of all pairs over the adjacency itself (n x n, diagonal 0).
    // Rows are done in bands: the l tiles of a band (the same 128 x 256
    // tiles, so the same values) go to a buffer first, then the band's upper
    // triangle is overwritten; later bands read only rows below it (rows of
    // finished bands that a tile touches give l values that are not used).
    // The lower triangle keeps the adjacency until it is mirrored at the end.
    // Extra memory: the band buffer, about n^2 / 16 values.
    void dissimilarity_inplace(bool similarity = false) {
        const int64_t N = n;
        const int64_t band_rt = std::max<int64_t>(1, N / 16 / kTileRows);
        const int n_rt = static_cast<int>((N + kTileRows - 1) / kTileRows);
        const int n_ct = static_cast<int>((N + kTileCols - 1) / kTileCols);
        std::vector<F> buf;
        for (int rt0 = 0; rt0 < n_rt; rt0 += static_cast<int>(band_rt)) {
            const int rt1 = std::min<int>(n_rt, rt0 + static_cast<int>(band_rt));
            const int ct0 = (rt0 * kTileRows) / kTileCols;
            const int64_t c_first = static_cast<int64_t>(ct0) * kTileCols;
            const int64_t W = N - c_first;  // buffer columns
            const int64_t r_first = static_cast<int64_t>(rt0) * kTileRows;
            const int64_t R_band = std::min<int64_t>(N, static_cast<int64_t>(rt1) * kTileRows)
                                   - r_first;
            buf.assign(static_cast<size_t>(R_band * W), F(0));
            struct Pair { int rt, ct; };
            std::vector<Pair> todo;
            for (int rt = rt0; rt < rt1; ++rt)
                for (int ct = (rt * kTileRows) / kTileCols; ct < n_ct; ++ct)
                    todo.push_back({rt, ct});
            parallel_for(pool, static_cast<int64_t>(todo.size()), [&](int64_t t) {
                const int r0 = todo[t].rt * kTileRows, R = std::min<int>(kTileRows, n - r0);
                const int c0 = todo[t].ct * kTileCols, C = std::min<int>(kTileCols, n - c0);
                const std::vector<F> L = tile(r0, c0);
                for (int a = 0; a < R; ++a)
                    std::copy(&L[static_cast<size_t>(a) * C], &L[static_cast<size_t>(a) * C] + C,
                              &buf[static_cast<size_t>((r0 + a - r_first) * W + (c0 - c_first))]);
            });
            parallel_for(pool, R_band, [&](int64_t a) {
                const int64_t i = r_first + a;
                const F *Lr = &buf[static_cast<size_t>(a * W)];
                F *row = &A[i * N];
                for (int64_t j = i + 1; j < N; ++j) {
                    const double l = Lr[j - c_first];
                    const double av = row[j];
                    const double t = (l + av) / (std::min(k[i], k[j]) + 1.0 - av);
                    row[j] = static_cast<F>(similarity ? t : 1.0 - t);
                }
                row[i] = similarity ? F(1) : F(0);
            });
        }
        parallel_for(pool, N, [&](int64_t j) {  // lower triangle <- upper
            for (int64_t i = 0; i < j; ++i) A[j * N + i] = A[i * N + j];
        });
    }
};

struct TomScorer::Impl : TomCore<float> {
    using TomCore<float>::TomCore;
};

TomScorer::TomScorer(const std::vector<std::vector<double>> &rows, int n_samples,
                     bool signed_net, double power, const std::vector<double> &powers,
                     WorkerPool &pool)
    : impl_(new Impl(pool)) {
    impl_->correlations(rows, n_samples);
    TomInfo info = impl_->choose_power(power, powers, n_samples, signed_net);
    power_ = info.power;
    estimated_ = info.estimated;
    fits_ = std::move(info.fits);
    impl_->adjacency(power_, signed_net);
}

TomScorer::~TomScorer() = default;

void TomScorer::score(const std::vector<CorEdge> &edges, std::vector<double> &tom) {
    Impl &m = *impl_;
    const size_t E = edges.size();
    tom.resize(E);
    const size_t n = static_cast<size_t>(m.n);
    for (size_t b = 0; b < E;) {
        const int rt = edges[b].i / kTileRows;
        size_t e = b;
        std::vector<int> cts;
        while (e < E && edges[e].i / kTileRows == rt) {
            const int c = edges[e].j / kTileCols;
            if (cts.empty() || cts.back() != c) {
                if (std::find(cts.begin(), cts.end(), c) == cts.end()) cts.push_back(c);
            }
            ++e;
        }
        if (rt < m.cur_row_tile) throw std::logic_error("TOM edges out of order.");
        m.ensure_tiles(rt, cts);
        const int r0 = rt * kTileRows;
        for (size_t x = b; x < e; ++x) {
            const int i = edges[x].i, j = edges[x].j;
            const int c = j / kTileCols, c0 = c * kTileCols;
            const int C = std::min(kTileCols, m.n - c0);
            const double l = m.tiles[c][static_cast<size_t>(i - r0) * C + (j - c0)];
            const double a = m.A[i * n + j];
            tom[x] = (l + a) / (std::min(m.k[i], m.k[j]) + 1.0 - a);
        }
        b = e;
    }
}

int TomScorer::n_genes() const { return impl_->n; }

TomInfo tom_dissimilarity(const std::vector<std::vector<double>> &rows, int n_samples,
                          bool signed_net, double power, const std::vector<double> &powers,
                          WorkerPool &pool, double *out, bool center, bool similarity) {
    TomCore<double> m(pool);
    m.correlations(rows, n_samples, out, center);
    TomInfo info = m.choose_power(power, powers, n_samples, signed_net);
    m.adjacency(info.power, signed_net);
    m.dissimilarity_inplace(similarity);
    return info;
}

std::vector<SoftThresholdFit> soft_threshold_fits(const std::vector<std::vector<double>> &rows,
                                                  int n_samples, bool signed_net,
                                                  const std::vector<double> &powers,
                                                  WorkerPool &pool, bool center) {
    TomCore<float> m(pool);  // as TomScorer picks the power
    m.correlations(rows, n_samples, nullptr, center);
    std::vector<double> pw(powers);
    std::sort(pw.begin(), pw.end());
    if (pw.empty() || pw.front() <= 0) throw std::invalid_argument("powers must be positive.");
    const auto kp = m.connectivity(pw, signed_net);
    std::vector<SoftThresholdFit> out;
    for (size_t p = 0; p < pw.size(); ++p) out.push_back(scale_free_fit(kp[p], pw[p]));
    return out;
}

namespace {

// WGCNA 1.74 pivot.c: pivot(v, len, target), v modified in place
double wgcna_pivot(double *v, size_t len, double target) {
    while (true) {
        if (len > 2) {
            size_t i1 = 0, i2 = len - 1, i3 = (len - 1) / 2, ip;
            if (v[i1] <= v[i2]) {
                if (v[i2] <= v[i3]) ip = i2;
                else if (v[i3] >= v[i1]) ip = i3;
                else ip = i1;
            } else {
                if (v[i1] <= v[i3]) ip = i1;
                else if (v[i2] <= v[i3]) ip = i3;
                else ip = i2;
            }
            double vp = v[ip];
            v[ip] = v[len - 1];
            v[len - 1] = vp;
            size_t bound = 0;
            for (size_t i = 0; i < len; i++)
                if (v[i] < vp) {
                    double x = v[bound];
                    v[bound] = v[i];
                    v[i] = x;
                    bound++;
                }
            v[len - 1] = v[bound];
            v[bound] = vp;
            double crit = target - bound;
            if (std::fabs(crit) > 1.0) {
                if (crit < 0) {
                    len = bound;
                } else {
                    v = v + bound + 1;
                    len = len - bound - 1;
                    target = target - bound - 1;
                }
                continue;
            }
            if (crit < 0) {
                double v1 = v[0];
                for (size_t i = 1; i < bound; i++) if (v[i] > v1) v1 = v[i];
                return v1 * (-crit) + vp * (1 + crit);
            }
            const double *w = v + bound + 1;
            const size_t lw = len - bound - 1;
            double v2 = w[0];
            for (size_t i = 1; i < lw; i++) if (w[i] < v2) v2 = w[i];
            return vp * (1 - crit) + v2 * crit;
        } else if (len == 2) {
            double v1 = std::min(v[0], v[1]), v2 = std::max(v[0], v[1]);
            if (target < 0) return v1;
            else if (target > 1) return v2;
            return target * v2 + (1 - target) * v1;
        }
        return v[0];
    }
}

}  // namespace

void wgcna_row_quantile(const std::vector<const double *> &sets, int64_t n, double q,
                        double *out, WorkerPool &pool) {
    const size_t ns = sets.size();
    const int64_t chunk = 1 << 16;
    const int64_t n_chunks = (n + chunk - 1) / chunk;
    parallel_for(pool, n_chunks, [&](int64_t c) {
        std::vector<double> v(ns);
        const int64_t a = c * chunk, b = std::min(n, a + chunk);
        for (int64_t i = a; i < b; ++i) {
            size_t m = 0;
            for (size_t s = 0; s < ns; ++s) {
                const double x = sets[s][i];
                if (!std::isnan(x)) v[m++] = x;  // NAs are left out, as WGCNA does
            }
            out[i] = m == 0 ? std::numeric_limits<double>::quiet_NaN()
                            : wgcna_pivot(v.data(), m, (1.0 * (m - 1)) * q);
        }
    });
}
