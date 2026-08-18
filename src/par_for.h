#ifndef STATMODELS7_PAR_FOR_H
#define STATMODELS7_PAR_FOR_H

#include <cstddef>
#include <RcppParallel.h>

// The parallel driver of the assembly kernels. The decomposition is over
// the elements of the OUTPUT: each element is accumulated in full by one
// thread, so no reduction is ever split and the result is bit-identical to
// the sequential order at any thread count (piano_parallel.txt, section 0,
// where the difference was measured at exactly zero). At threads <= 1
// nothing parallel is entered; the thresholds below which a caller stays
// sequential are internal and measured, not arguments.
namespace sm7 {

template <typename Body>
struct BodyWorker : public RcppParallel::Worker {
  const Body& body;
  explicit BodyWorker(const Body& b) : body(b) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) body(i);
  }
};

template <typename Body>
inline void par_for(std::size_t n, int threads, const Body& body) {
  if (threads > 1) {
    BodyWorker<Body> w(body);
    RcppParallel::parallelFor(0, n, w);
  } else {
    for (std::size_t i = 0; i < n; ++i) body(i);
  }
}

} // namespace sm7

#endif
