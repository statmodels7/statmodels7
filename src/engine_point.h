// The point of the outer search (phase 4): what evaluate() computes at one
// outer vector, behind an interface that engine_search.cpp (the search, on
// Armadillo) and engine_iwls.cpp (the compiled point, on the engine's own
// matrices) share without either seeing the other's headers.
#ifndef STATMODELS7_ENGINE_POINT_H
#define STATMODELS7_ENGINE_POINT_H

#include <Rcpp.h>
#include <string>
#include <vector>

// what the point source returns for one evaluation
struct PointOut {
    bool ok = false;
    std::string err;       // the message of an error at a later evaluation
    bool has_err = false;
    double value = NA_REAL, loglik = NA_REAL, penalty = NA_REAL;
    std::vector<double> par;   // the mode
    std::vector<double> vals;  // v_value() per label
    std::string hv;            // the unavailable line's values
    bool hv_num = false;       // the line is formatted in R from `vals`
    std::string cause;         // the first evaluation's message, when unavailable
};

// R's evaluate() body, criterion_resolution() and the derivatives at the
// last point evaluated
class PointSource {
public:
    virtual ~PointSource() {}
    virtual void eval(const std::vector<double>& v, const std::vector<double>& warm,
                      bool first, PointOut& out) = 0;
    // false where the derivative is NULL
    virtual bool grad(std::vector<double>& g) = 0;
    virtual bool hess(std::vector<double>& H) = 0;
    // criterion_resolution() at the last point: r (NA if none) and, where it
    // was refused on the mode error, that error
    virtual void reading(double& r, bool& has_me, double& me) = 0;
    // how many points the compiled source handed to `fallback`
    virtual int handed_back() const { return 0; }
    // the final refit at v from warm (D40): the inner fit's raw result, the
    // objective before and after, and the criterion; an empty list where the
    // refit runs in R
    virtual Rcpp::List refit(const std::vector<double>& v, const std::vector<double>& warm) {
        (void) v; (void) warm;
        return Rcpp::List();
    }
};

// the compiled point (engine_iwls.cpp); `fallback` computes a point that the
// compiled one hands back, and is owned by the caller
PointSource* make_compiled_point_source(Rcpp::List cfg, Rcpp::List cb, PointSource* fallback);

#endif
