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
- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  : Estimate the Hyperparameters by a Marginal Likelihood
- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  : Prediction-Error Criteria for the Hyperparameters
- [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) :
  Choose the Hyperparameters by Cross-Validation
- [`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
  : Simulate a Response From a Written Model

## Reading a fit

The model as a function of parameters and data, any parameter or moment
predicted at new data, what the fit says about its own uncertainty, and
the elapsed time in the unit it deserves.

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
- [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  : How the Hyperparameters Are Estimated
- [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  : A Fitted Model
- [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  : The Specification of a Model, Before It Is Fitted
- [`StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md)
  : A Summary of a Fitted Model
- [`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
  : Solve a Scoring Step From the Square-Root Design
- [`block_predictors()`](https://statmodels7.github.io/statmodels7/reference/block_predictors.md)
  : The Predictors a Coefficient Direction Induces
- [`blocks_at_kink()`](https://statmodels7.github.io/statmodels7/reference/blocks_at_kink.md)
  : Record Where a Path Has Just Been
- [`bounded_bump()`](https://statmodels7.github.io/statmodels7/reference/bounded_bump.md)
  : Move a Hyperparameter Without Leaving Its Interval
- [`check_offsets()`](https://statmodels7.github.io/statmodels7/reference/check_offsets.md)
  : Validate Offsets
- [`check_weights()`](https://statmodels7.github.io/statmodels7/reference/check_weights.md)
  : Validate Prior Weights
- [`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
  : Cholesky Factors of Small Blocks, Vectorized Over Observations
- [`coef(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
  : The Coefficients of a Fit
- [`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
  : Where Each Stacked Coefficient Comes From
- [`confint(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  : Confidence Intervals for a Fit
- [`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
  : The Third Derivative of the Objective Contracted Once
- [`contract4()`](https://statmodels7.github.io/statmodels7/reference/contract4.md)
  : The Fourth Derivative of the Objective Contracted Twice
- [`coord_covariance()`](https://statmodels7.github.io/statmodels7/reference/coord_covariance.md)
  : Which Way of Holding the Gradient Is Cheaper
- [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
  : Fit a Separable Block by Coordinate Descent
- [`coord_offset()`](https://statmodels7.github.io/statmodels7/reference/coord_offset.md)
  : The Offset of One Equation
- [`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md)
  : Which Coordinates a Path Point Has to Visit
- [`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
  : The Working Response and Weights of One Equation
- [`criterion_tol()`](https://statmodels7.github.io/statmodels7/reference/criterion_tol.md)
  : The Tolerance a Criterion Asks For
- [`cv_curve()`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md)
  : The Held-Out Deviance of Every Point of a Path
- [`cv_folds()`](https://statmodels7.github.io/statmodels7/reference/cv_folds.md)
  : Assign the Observations to Folds
- [`d3_key()`](https://statmodels7.github.io/statmodels7/reference/d3_key.md)
  : The Name of a Third-Derivative Component
- [`d4_key()`](https://statmodels7.github.io/statmodels7/reference/d4_key.md)
  : The Name of a Fourth-Derivative Component
- [`draw_coefficients()`](https://statmodels7.github.io/statmodels7/reference/draw_coefficients.md)
  : Draw or Validate the Coefficients of a Simulation
- [`fit_expected()`](https://statmodels7.github.io/statmodels7/reference/fit_expected.md)
  : Which Information Matrix a Fit Used
- [`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
  : Fit the Smooth Block
- [`fitted(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
  : The Fitted Values of a Model
- [`fitted_ranges()`](https://statmodels7.github.io/statmodels7/reference/fitted_ranges.md)
  : What Each Distribution Parameter Reached
- [`flat_directions()`](https://statmodels7.github.io/statmodels7/reference/flat_directions.md)
  : Which Coefficients a Singular Curvature Is Flat In
- [`fmt_step()`](https://statmodels7.github.io/statmodels7/reference/fmt_step.md)
  : Format a Step Length
- [`format_conflicts()`](https://statmodels7.github.io/statmodels7/reference/format_conflicts.md)
  : Render a Conflict Report
- [`free_scale()`](https://statmodels7.github.io/statmodels7/reference/free_scale.md)
  [`free_scale2()`](https://statmodels7.github.io/statmodels7/reference/free_scale.md)
  [`link_slopes()`](https://statmodels7.github.io/statmodels7/reference/free_scale.md)
  : Carry a Hyperparameter Derivative onto the Free Scale
- [`hess_key()`](https://statmodels7.github.io/statmodels7/reference/hess_key.md)
  : The Name of a Second-Derivative Component
- [`hyper_key()`](https://statmodels7.github.io/statmodels7/reference/hyper_key.md)
  : Resolve a Term's Name Against a Specification
- [`hyper_plain()`](https://statmodels7.github.io/statmodels7/reference/hyper_plain.md)
  : The Hyperparameters in the Shape statmod() Accepts
- [`hyper_set()`](https://statmodels7.github.io/statmodels7/reference/hyper_set.md)
  : Set One Hyperparameter
- [`hyper_to_eta()`](https://statmodels7.github.io/statmodels7/reference/hyper_to_eta.md)
  [`eta_to_hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper_to_eta.md)
  : Move Between the Hyperparameters and the Free Vector
- [`hyper_value()`](https://statmodels7.github.io/statmodels7/reference/hyper_value.md)
  : One Hyperparameter of an Index
- [`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
  : The Per-Observation Information Blocks
- [`inner_settings()`](https://statmodels7.github.io/statmodels7/reference/inner_settings.md)
  : What the Inner Method Says About How to Fit
- [`integrated_basis()`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md)
  : The Subspace a Marginal Criterion Integrates Over
- [`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
  : Fit the Smooth Block by Iterated Weighted Least Squares
- [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  : The Pieces One Scoring Step Needs
- [`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  : Solve One Weighted Least Squares Step
- [`kink_hypers()`](https://statmodels7.github.io/statmodels7/reference/kink_hypers.md)
  : Which Hyperparameters Set the Size of the Kink
- [`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)
  : The Size of a Penalty's Kink
- [`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
  : The Hyperparameter That Gives the Kink a Chosen Size
- [`logLik(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md)
  : The Maximized Log-Likelihood of a Fit
- [`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`gradient()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  [`hessian()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
  : The Model as a Function of Parameters and Data
- [`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
  : The Budget and the Stopping Rule of the Alternation
- [`one_sided()`](https://statmodels7.github.io/statmodels7/reference/one_sided.md)
  : Build a One-Sided Formula From an Expression
- [`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
  : Estimate the Hyperparameters
- [`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
  : Can the Exact Gradient Be Computed Here?
- [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
  : The Hyperparameters an Outer Method Estimates
- [`outer_k()`](https://statmodels7.github.io/statmodels7/reference/outer_k.md)
  : The Price of One Degree of Freedom
- [`outer_minimize()`](https://statmodels7.github.io/statmodels7/reference/outer_minimize.md)
  : Is a Criterion Minimized?
- [`outer_path_defaults()`](https://statmodels7.github.io/statmodels7/reference/outer_path_defaults.md)
  : The Defaults a Path Carries
- [`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
  : The Per-Hyperparameter Pieces of the Outer Derivatives
- [`outer_tau()`](https://statmodels7.github.io/statmodels7/reference/outer_tau.md)
  : The Effective Degrees of Freedom of a Whole Fit
- [`pair_key()`](https://statmodels7.github.io/statmodels7/reference/pair_key.md)
  : The Key of a Hyperparameter Pair
- [`par_at()`](https://statmodels7.github.io/statmodels7/reference/par_at.md)
  : Resolve a Parameter Structure
- [`path_block()`](https://statmodels7.github.io/statmodels7/reference/path_block.md)
  : The Block a Path Row Belongs To
- [`path_null_score()`](https://statmodels7.github.io/statmodels7/reference/path_null_score.md)
  : The Largest Score a Kinked Block Has to Beat
- [`path_pick()`](https://statmodels7.github.io/statmodels7/reference/path_pick.md)
  : Choose a Point of the Path
- [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
  : Which Hyperparameters a Path Has to Select
- [`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
  : The Values a Path Visits
- [`pd_repair()`](https://statmodels7.github.io/statmodels7/reference/pd_repair.md)
  : Floor the Eigenvalues of a Curvature Matrix
- [`penalty_answers()`](https://statmodels7.github.io/statmodels7/reference/penalty_answers.md)
  : Does a Penalty Supply What a Marginal Criterion Needs?
- [`penalty_has_kink()`](https://statmodels7.github.io/statmodels7/reference/penalty_has_kink.md)
  : Does a Penalty Have a Kink?
- [`penalty_range_basis()`](https://statmodels7.github.io/statmodels7/reference/penalty_range_basis.md)
  : An Orthonormal Basis of a Penalty's Range Space
- [`penalty_sqrt()`](https://statmodels7.github.io/statmodels7/reference/penalty_sqrt.md)
  : A Square-Root Factor of the Penalty
- [`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
  : A Penalty's Starting Hyperparameters
- [`predict(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
  : Predict From a Fitted Model
- [`predict_moments()`](https://statmodels7.github.io/statmodels7/reference/predict_moments.md)
  : The Quantities a Fit Can Predict
- [`print(`*`<OuterMethod>`*`)`](https://statmodels7.github.io/statmodels7/reference/print.OuterMethod.md)
  : Print an Outer Method
- [`print(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md)
  : Print a Fitted Model
- [`print(`*`<StatmodSummary>`*`)`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSummary.md)
  : Print a Model Summary
- [`print_block()`](https://statmodels7.github.io/statmodels7/reference/print_block.md)
  : Print One Block of a Summary
- [`refreshes_own_block()`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md)
  : Does a Term Recompute Its Own Block?
- [`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md)
  : Reject a Term the Fitting Scheme Does Not Cover
- [`smooth_linear_cols()`](https://statmodels7.github.io/statmodels7/reference/smooth_linear_cols.md)
  : Which Coefficients of a Smooth Are the Linear Part
- [`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  : Invert a Matrix That Ought to Be Positive Definite
- [`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
  : Fit One Non-Smooth Block, the Others Held Fixed
- [`spec_at()`](https://statmodels7.github.io/statmodels7/reference/spec_at.md)
  : Rebuild a Specification Against New Data
- [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  : The Square-Root Design
- [`statmod_active()`](https://statmodels7.github.io/statmodels7/reference/statmod_active.md)
  : Which Coefficients Are Not Sitting at a Kink
- [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  : The Alternation Between the Smooth Block and the Rest
- [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  : Split a Specification's Terms Into the Smooth Block and the Rest
- [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
  : Effective Degrees of Freedom, Per Term
- [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  : The Linear Predictors and the Parameters They Give
- [`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md)
  : Override the Starting Hyperparameters
- [`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
  : The Hyperparameters a Specification Starts From
- [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
  : The Information of the Weighted Log-Likelihood
- [`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
  : The Intercept of Each Equation, on the Link Scale
- [`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
  : The Weighted Log-Likelihood of a Specification at Given Coefficients
- [`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
  : The Marginal Criterion at Given Coefficients and Hyperparameters
- [`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
  : The Exact Gradient of the Marginal Criterion
- [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
  : The Exact Hessian of the Marginal Criterion
- [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  : The Objective, Its Gradient and Its Hessian, Stacked
- [`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
  : Select the Hyperparameters of a Kinked Penalty Along a Path
- [`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)
  : A Prediction-Error Criterion at Given Coefficients and
  Hyperparameters
- [`statmod_pe_derivs()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe_derivs.md)
  : The Exact Derivatives of a Prediction-Error Criterion
- [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  : Every Penalized Unit of a Specification
- [`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
  : The Penalty of a Specification at Given Coefficients
- [`statmod_penalty_keys()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_keys.md)
  : Every Penalty in a Model, Without the Design
- [`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
  : The Same Model Read on Other Rows
- [`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
  : The Score of the Weighted Log-Likelihood
- [`statmod_select()`](https://statmodels7.github.io/statmodels7/reference/statmod_select.md)
  : Estimate the Hyperparameters, by Whichever Route Each One Admits
- [`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
  : Starting Coefficients
- [`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
  : Interpret and Build Each Parameter's Terms
- [`statmod_unit()`](https://statmodels7.github.io/statmodels7/reference/statmod_unit.md)
  : One Penalized Unit, by Parameter and Key
- [`statmodels7`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  [`statmodels7-package`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  : statmodels7: The S7 Toolkit for Statistical Modeling
- [`statmodels7_attach()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach.md)
  : Attach the Member Packages
- [`statmodels7_attach_message()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach_message.md)
  : The Body of the Attach Message
- [`summary(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  : Summarize a Fitted Model
- [`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
  : The Blocks of One Distribution Parameter
- [`term_block_kind()`](https://statmodels7.github.io/statmodels7/reference/term_block_kind.md)
  : What Kind of Block a Term Reports As
- [`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
  : Evaluate a Formula's Terms With modelterms7 in Front
- [`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
  : The Trace of the Determinant's Movement With the Coefficients
- [`unfittable_reason()`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md)
  : Why a Term Is Outside the Fitting Scheme
- [`unknown_what()`](https://statmodels7.github.io/statmodels7/reference/unknown_what.md)
  : The Message for an Unrecognized Prediction Target
- [`vb_inner()`](https://statmodels7.github.io/statmodels7/reference/vb_inner.md)
  : The Verbosity of an Inner Fit Inside the Outer Search
- [`vcov(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  : The Variance Matrix of a Fit
- [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  : Resolve the Verbosity Setting
