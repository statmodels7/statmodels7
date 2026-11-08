// joint_decrement() in a translation unit of its own: it calls R's LAPACK
// through R_ext/Lapack.h, whose prototypes clang rejects beside
// RcppArmadillo's in one translation unit (macOS).
#ifndef SM7_ENGINE_SEARCH_LA_H
#define SM7_ENGINE_SEARCH_LA_H

#include <vector>

double joint_decrement(const std::vector<double>& g, const std::vector<double>& A);

#endif
