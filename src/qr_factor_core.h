// The Householder factor of qr_factor_cpp(), on a column-major m x p buffer
// overwritten in place; R (p x p, upper) is written to `out`. Shared by
// qr_factor_cpp() and the compiled inner fit, so that both routes run the same
// arithmetic.
#ifndef SM7_QR_FACTOR_CORE_H
#define SM7_QR_FACTOR_CORE_H

#include <R.h>
#include <vector>
#include <cmath>
#include "par_for.h"

namespace sm7 {

inline void qr_factor_core(double* Bp, int m, int p, int threads, double* out) {
    std::vector<double> v(m > 0 ? m : 1);
    const int kmax = (m < p) ? m : p;

    // measured on this machine: the per-step region costs ~44 us to open,
    // and a trailing update runs at roughly one element per nanosecond, so
    // a step is worth splitting from about a quarter of a million elements
    const double kMinWork = 2.6e5;

    for (int j = 0; j < kmax; ++j) {
        const int len = m - j;
        double* col = Bp + (std::size_t) j * m + j;
        double nrm = 0.0;
        for (int i = 0; i < len; ++i) nrm += col[i] * col[i];
        nrm = std::sqrt(nrm);
        if (!R_FINITE(nrm) || nrm == 0.0) continue;
        const double alpha = (col[0] >= 0.0) ? -nrm : nrm;
        for (int i = 0; i < len; ++i) v[i] = col[i];
        v[0] -= alpha;
        double vtv = 0.0;
        for (int i = 0; i < len; ++i) vtv += v[i] * v[i];
        if (vtv == 0.0) continue;
        col[0] = alpha;
        for (int i = 1; i < len; ++i) col[i] = 0.0;

        const int ntrail = p - j - 1;
        if (ntrail <= 0) continue;
        const double* vp = v.data();
        const double inv = 2.0 / vtv;
        auto body = [&](std::size_t t) {
            const int k = j + 1 + (int) t;
            double* ck = Bp + (std::size_t) k * m + j;
            double s = 0.0;
            for (int i = 0; i < len; ++i) s += vp[i] * ck[i];
            s *= inv;
            for (int i = 0; i < len; ++i) ck[i] -= s * vp[i];
        };
        const double work = (double) ntrail * (double) len;
        sm7::par_for((std::size_t) ntrail, (work >= kMinWork) ? threads : 1,
                     body);
    }

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) out[(std::size_t) j * p + i] = 0.0;
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i <= j && i < kmax; ++i) {
            out[(std::size_t) j * p + i] = Bp[(std::size_t) j * m + i];
        }
    }
}

}  // namespace sm7

#endif
