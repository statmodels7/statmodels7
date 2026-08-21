#include <Rcpp.h>
#include "par_for.h"

using namespace Rcpp;

// The per-observation leverage diagonal of one block pair, over the nonzeros
// of two sparse designs:
//
//   G_i = sum_{p in row i of A} sum_{q in row i of B} A_ip B_iq M[j_p, k_q]
//
// Where a design is built from grouping indicators a row has a handful of
// nonzeros, so this is a short double loop per observation -- a scalar
// recursion, which the R form could only express as seven vectors as long as
// the PAIR count: at 3.38 million pairs that is about 190 MB allocated per
// call, and the call was measured at 262 ms and 31.6 per cent of a REML fit
// over 500 random-effect levels.
//
// The decomposition is over the elements of the OUTPUT: row i's sum is
// computed and written in full by one thread, in the order the sequential
// loop performs it, so the result does not depend on the thread count. The
// body reads and writes preallocated memory and calls nothing.
//
// Both designs arrive ordered by row, with the start and the count of each
// row's nonzeros, so a row's two ranges are contiguous.
// [[Rcpp::export]]
NumericVector leverage_pairs_cpp(IntegerVector astart, IntegerVector acnt,
                                 IntegerVector aj, NumericVector av,
                                 IntegerVector bstart, IntegerVector bcnt,
                                 IntegerVector bj, NumericVector bv,
                                 NumericMatrix M, int n, int threads) {
    NumericVector out(n);
    const int* asp = astart.begin();
    const int* acp = acnt.begin();
    const int* ajp = aj.begin();
    const double* avp = av.begin();
    const int* bsp = bstart.begin();
    const int* bcp = bcnt.begin();
    const int* bjp = bj.begin();
    const double* bvp = bv.begin();
    const double* Mp = M.begin();
    const int mrows = M.nrow();
    double* op = out.begin();

    // one region is worth opening from a few thousand rows: the body is a
    // handful of multiply-adds an observation, the cheap class of par_for.h
    sm7::par_for((std::size_t) n, (n >= 4096) ? threads : 1,
                 [&](std::size_t i) {
        const int na = acp[i], nb = bcp[i];
        if (na == 0 || nb == 0) return;
        const int a0 = asp[i], b0 = bsp[i];
        double s = 0.0;
        for (int p = 0; p < na; ++p) {
            const double va = avp[a0 + p];
            const int jr = ajp[a0 + p] - 1;          // M is column-major
            for (int q = 0; q < nb; ++q) {
                s += va * bvp[b0 + q] *
                     Mp[jr + (std::size_t) (bjp[b0 + q] - 1) * mrows];
            }
        }
        op[i] = s;
    });
    return out;
}
