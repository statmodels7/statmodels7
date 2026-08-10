# Package index

## The toolkit

Which packages the toolkit is made of, at which versions, and whether
any of their exports mask one another.

- [`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
  : The Packages of the Toolkit
- [`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
  : Installed Versions of the Toolkit
- [`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
  : Exports of the Toolkit That Mask One Another

## Specifying a model

One formula carries every parameter of the distribution, the equations
separated by a bar, and each equation’s terms come from modelterms7.

- [`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)
  : Build a Model Specification
- [`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
  : Split a Multi-Parameter Formula Into One Equation Per Parameter
- [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  : The Design of a Specification

## Fitting

One formula, one distribution, and a scheme that fits the terms whose
penalties are differentiable all together while the others are
alternated around them.

- [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  : Fit a Model
- [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  [`print(`*`<Iwls>`*`)`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  : Iterated Weighted Least Squares
- [`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
  : Simulate a Response From a Written Model

## Reading a fit

The model as a function of parameters and data, any parameter or moment
predicted at new data, and the elapsed time in the unit it deserves.

- [`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`gradient()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`hessian()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  : The Model as a Function of Parameters and Data
- [`format_duration()`](https://statmodels7.github.io/statmodels7/reference/format_duration.md)
  : Format a Duration in the Unit It Deserves

## Installing and updating

The toolkit is not on CRAN, so the install path is GitHub and the
dependencies among the members are resolved by pak.

- [`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
  : Install or Update the Toolkit

## Internals

The machinery the exported functions are built from. None of it is
exported, and none is needed to use the package.

- [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  : The Iterated Weighted Least Squares Method
- [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  : A Fitted Model
- [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  : The Specification of a Model, Before It Is Fitted
- [`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
  : Solve a Scoring Step From the Square-Root Design
- [`check_offsets()`](https://statmodels7.github.io/statmodels7/reference/check_offsets.md)
  : Validate Offsets
- [`check_weights()`](https://statmodels7.github.io/statmodels7/reference/check_weights.md)
  : Validate Prior Weights
- [`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
  : Cholesky Factors of Small Blocks, Vectorized Over Observations
- [`coef(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
  : The Coefficients of a Fit
- [`draw_coefficients()`](https://statmodels7.github.io/statmodels7/reference/draw_coefficients.md)
  : Draw or Validate the Coefficients of a Simulation
- [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  : Fit the Smooth Block
- [`fitted(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
  : The Fitted Values of a Model
- [`fmt_step()`](https://statmodels7.github.io/statmodels7/reference/fmt_step.md)
  : Format a Step Length
- [`format_conflicts()`](https://statmodels7.github.io/statmodels7/reference/format_conflicts.md)
  : Render a Conflict Report
- [`hess_key()`](https://statmodels7.github.io/statmodels7/reference/hess_key.md)
  : The Name of a Second-Derivative Component
- [`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
  : The Per-Observation Information Blocks
- [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
  : Fit the Smooth Block by Iterated Weighted Least Squares
- [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  : The Pieces One Scoring Step Needs
- [`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  : Solve One Weighted Least Squares Step
- [`logLik(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
  : The Maximized Log-Likelihood of a Fit
- [`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`gradient()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`hessian()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  : The Model as a Function of Parameters and Data
- [`one_sided()`](https://statmodels7.github.io/statmodels7/reference/one_sided.md)
  : Build a One-Sided Formula From an Expression
- [`par_at()`](https://statmodels7.github.io/statmodels7/reference/par_at.md)
  : Resolve a Parameter Structure
- [`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
  : Floor the Eigenvalues of a Curvature Matrix
- [`penalty_has_kink()`](https://statmodels7.github.io/statmodels7/reference/penalty_has_kink.md)
  : Does a Penalty Have a Kink?
- [`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
  : A Square-Root Factor of the Penalty
- [`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
  : A Penalty's Starting Hyperparameters
- [`predict(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
  : Predict From a Fitted Model
- [`predict_moments()`](https://statmodels7.github.io/statmodels7/reference/predict_moments.md)
  : The Quantities a Fit Can Predict
- [`print(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md)
  : Print a Fitted Model
- [`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
  : Fit One Non-Smooth Block, the Others Held Fixed
- [`spec_at()`](https://statmodels7.github.io/statmodels7/reference/spec_at.md)
  : Rebuild a Specification Against New Data
- [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  : The Square-Root Design
- [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  : Split a Specification's Terms Into the Smooth Block and the Rest
- [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
  : Effective Degrees of Freedom, Per Term
- [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  : The Linear Predictors and the Parameters They Give
- [`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
  : The Hyperparameters a Specification Starts From
- [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
  : The Information of the Weighted Log-Likelihood
- [`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
  : The Weighted Log-Likelihood of a Specification at Given Coefficients
- [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  : The Objective, Its Gradient and Its Hessian, Stacked
- [`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
  : The Penalty of a Specification at Given Coefficients
- [`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
  : The Score of the Weighted Log-Likelihood
- [`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
  : Starting Coefficients
- [`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
  : Interpret and Build Each Parameter's Terms
- [`statmodels7`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  [`statmodels7-package`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  : statmodels7: The S7 Toolkit for Statistical Modeling
- [`statmodels7_attach()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach.md)
  : Attach the Member Packages
- [`statmodels7_attach_message()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach_message.md)
  : The Body of the Attach Message
- [`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
  : Evaluate a Formula's Terms With modelterms7 in Front
- [`unknown_what()`](https://statmodels7.github.io/statmodels7/reference/unknown_what.md)
  : The Message for an Unrecognized Prediction Target
- [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  : Resolve the Verbosity Setting
