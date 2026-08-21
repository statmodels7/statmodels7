#ifndef STATMODELS7_PAR_FOR_H
#define STATMODELS7_PAR_FOR_H

#include <cstddef>
#include <fenv.h>
#include <RcppParallel.h>

// The parallel driver of the assembly kernels. The decomposition is over
// the elements of the OUTPUT: each element is accumulated in full by one
// thread, so no reduction is ever split and the result is bit-identical to
// the sequential order at any thread count (piano_parallel.txt, section 0,
// where the difference was measured at exactly zero). At threads <= 1
// nothing parallel is entered; the thresholds below which a caller stays
// sequential are internal and measured, not arguments.
//
// The shape is distributions7's d7::par_for(), which arrived at it the
// expensive way, and this header did not follow until 2026-08-21. Two
// things it now shares:
//
//   - the loop is NOINLINE and the sequential branch runs THROUGH the
//     worker rather than through a loop of its own, so the two branches
//     execute one compiled copy. An inlined sequential loop and a worker's
//     out-of-line call are two bodies the compiler is free to optimize
//     apart, which is what produced one-ulp cross-count differences in
//     distributions7 on the Windows CI runner while nothing about the
//     decomposition had changed. These kernels are dot products of plain
//     doubles and the identity has held here, but it was holding by the
//     optimizer's leave rather than by construction;
//
//   - the worker installs the calling thread's floating-point environment
//     before its chunk. Nothing here calls an Rmath routine, so this buys
//     no correctness today; it costs one call per chunk and makes the rule
//     uniform across the toolkit's three drivers, so a kernel added later
//     cannot acquire the defect by being written in the obvious way.
namespace sm7 {

#if defined(__GNUC__) || defined(__clang__)
#define SM7_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define SM7_NOINLINE __declspec(noinline)
#else
#define SM7_NOINLINE
#endif

template <typename Body>
struct BodyWorker : public RcppParallel::Worker {
  const Body& body;
  fenv_t env;
  explicit BodyWorker(const Body& b) : body(b) { fegetenv(&env); }
  SM7_NOINLINE void operator()(std::size_t begin, std::size_t end) {
    fesetenv(&env);
    for (std::size_t i = begin; i < end; ++i) body(i);
  }
};

// The count is passed to parallelFor rather than left to the process-level
// setting: resolveValue() prefers an explicit positive value to
// RCPP_PARALLEL_NUM_THREADS, so a fit that sized the pool through
// numericals7::local_threads() keeps that size and a caller that did not
// gets the count it asked for instead of every core the machine has.
template <typename Body>
inline void par_for(std::size_t n, int threads, const Body& body) {
  BodyWorker<Body> w(body);
  if (threads > 1) {
    RcppParallel::parallelFor(0, n, w, 1, threads);
  } else {
    w(0, n);
  }
}

} // namespace sm7

#endif
