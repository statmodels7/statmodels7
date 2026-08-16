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
- [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  : Options for the Unpenalized Parametric Block
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

## Where a fit begins

A starting value as a strategy rather than a vector of numbers, asked
once before the fit alternates. The default is the intercept-only fit;
what the others offer is a procedure, from a random perturbation to a
global search of the likelihood.

- [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
  : S7 Class for a Starting-Value Strategy
- [`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
  : Where a Fit Begins
- [`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
  : Start at the Intercept-Only Fit
- [`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md)
  : Start at the Origin
- [`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
  : Start From a Random Draw
- [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
  : Search the Likelihood for a Starting Point

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
- [`StartIntercepts()`](https://statmodels7.github.io/statmodels7/reference/StartIntercepts-class.md)
  [`StartOrigin()`](https://statmodels7.github.io/statmodels7/reference/StartIntercepts-class.md)
  [`StartRandom()`](https://statmodels7.github.io/statmodels7/reference/StartIntercepts-class.md)
  [`StartSearch()`](https://statmodels7.github.io/statmodels7/reference/StartIntercepts-class.md)
  : S7 Classes for the Shipped Strategies
- [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  : A Fitted Model
- [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  : The Specification of a Model, Before It Is Fitted
- [`StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md)
  : A Summary of a Fitted Model
- [`add_offsets()`](https://statmodels7.github.io/statmodels7/reference/add_offsets.md)
  : Add Two Sets of Offsets
- [`answers_term_third()`](https://statmodels7.github.io/statmodels7/reference/answers_term_third.md)
  : Does a Term Supply Its Third Derivative?
- [`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md)
  : Solve a Scoring Step From the Square-Root Design
- [`bind_blocks()`](https://statmodels7.github.io/statmodels7/reference/bind_blocks.md)
  : Bind a Model's Term Blocks Side by Side
- [`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
  : The Per-Observation Diagonal of Each Block of a Matrix
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
- [`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
  : The Penalized Block, in the Storage It Arrived In
- [`coord_call()`](https://statmodels7.github.io/statmodels7/reference/coord_call.md)
  : Run the Compiled Coordinate Descent on Either Storage
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
- [`ctx_deriv()`](https://statmodels7.github.io/statmodels7/reference/ctx_deriv.md)
  : A Higher Derivative of the Log-Density at the Context's Point
- [`ctx_information()`](https://statmodels7.github.io/statmodels7/reference/ctx_information.md)
  : The Information at the Context's Point
- [`ctx_leverage()`](https://statmodels7.github.io/statmodels7/reference/ctx_leverage.md)
  : The Per-Observation Diagonals at the Context's Point
- [`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
  : The Penalized Matrix and Its Inverse
- [`ctx_penalty()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalty.md)
  : The Penalty's Hessian at the Context's Point
- [`ctx_theta()`](https://statmodels7.github.io/statmodels7/reference/ctx_theta.md)
  : The Linear Predictors' Parameters at the Context's Point
- [`ctx_trace_matrix()`](https://statmodels7.github.io/statmodels7/reference/ctx_trace_matrix.md)
  : The Matrix the Traces Are Taken Against
- [`ctx_usable()`](https://statmodels7.github.io/statmodels7/reference/ctx_usable.md)
  : Refuse a Context That Belongs Somewhere Else
- [`cv_bind_inputs()`](https://statmodels7.github.io/statmodels7/reference/cv_bind_inputs.md)
  : Carry a Term's Matrix Input Onto a Subset of the Rows
- [`cv_curve()`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md)
  : The Held-Out Deviance of Every Point of a Path
- [`cv_folds()`](https://statmodels7.github.io/statmodels7/reference/cv_folds.md)
  : Assign the Observations to Folds
- [`d3_key()`](https://statmodels7.github.io/statmodels7/reference/d3_key.md)
  : The Name of a Third-Derivative Component
- [`d4_key()`](https://statmodels7.github.io/statmodels7/reference/d4_key.md)
  : The Name of a Fourth-Derivative Component
- [`deriv3_key()`](https://statmodels7.github.io/statmodels7/reference/deriv3_key.md)
  : The Name of a Third-Derivative Component
- [`deriv4_key()`](https://statmodels7.github.io/statmodels7/reference/deriv4_key.md)
  : The Name of a Fourth-Derivative Component
- [`design_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  [`as_dense()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  [`as_sparse()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  [`zero_information()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md)
  : Is a Design Sparse, and the Zero Matrix to Accumulate It Into
- [`diagonal_sqrt()`](https://statmodels7.github.io/statmodels7/reference/diagonal_sqrt.md)
  : The Factor of a Diagonal Penalty
- [`.structural_blocks()`](https://statmodels7.github.io/statmodels7/reference/dot-structural_blocks.md)
  : The Model's Derivative Pieces for a Filter's Recursion
- [`draw_coefficients()`](https://statmodels7.github.io/statmodels7/reference/draw_coefficients.md)
  : Draw or Validate the Coefficients of a Simulation
- [`eval_offsets()`](https://statmodels7.github.io/statmodels7/reference/eval_offsets.md)
  : Evaluate the Offsets a Formula Names
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
- [`iwls_info_diag()`](https://statmodels7.github.io/statmodels7/reference/iwls_info_diag.md)
  : The Diagonal of the Information a Step Uses
- [`iwls_met()`](https://statmodels7.github.io/statmodels7/reference/iwls_met.md)
  : Has the Step's Stopping Rule Been Met?
- [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  : The Pieces One Scoring Step Needs
- [`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
  : The Dimensionless Reading of the Stopping Rule
- [`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  : Solve One Weighted Least Squares Step
- [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md)
  : The Rows of the Joint Predictor Derivative
- [`kink_by_power()`](https://statmodels7.github.io/statmodels7/reference/kink_by_power.md)
  : Invert the Size of the Kink Through a Power Law
- [`kink_hypers()`](https://statmodels7.github.io/statmodels7/reference/kink_hypers.md)
  : Which Hyperparameters Set the Size of the Kink
- [`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md)
  : How the Size of the Kink Scales With a Hyperparameter
- [`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)
  : The Size of a Penalty's Kink
- [`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
  : The Hyperparameter That Gives the Kink a Chosen Size
- [`leverage_pairs()`](https://statmodels7.github.io/statmodels7/reference/leverage_pairs.md)
  : The Leverage Diagonal Over the Nonzeros of Two Rows
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
- [`outer_context()`](https://statmodels7.github.io/statmodels7/reference/outer_context.md)
  : One Evaluation Point, Shared
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
  : The Properties Every Criterion Carries
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
- [`path_bounded()`](https://statmodels7.github.io/statmodels7/reference/path_bounded.md)
  : Is a Hyperparameter Bounded Above?
- [`path_by_kink()`](https://statmodels7.github.io/statmodels7/reference/path_by_kink.md)
  : Is a Hyperparameter Swept by the Size of Its Kink?
- [`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
  : What a Path Does Where the Term Says Nothing
- [`path_forced()`](https://statmodels7.github.io/statmodels7/reference/path_forced.md)
  : The Values a Caller Wrote Out
- [`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md)
  : The Values a Path Visits Over a Bounded Hyperparameter
- [`path_null_score()`](https://statmodels7.github.io/statmodels7/reference/path_null_score.md)
  : The Largest Score a Kinked Block Has to Beat
- [`path_pick()`](https://statmodels7.github.io/statmodels7/reference/path_pick.md)
  : Choose a Point of the Path
- [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)
  : Which Hyperparameters a Path Has to Select
- [`path_steps()`](https://statmodels7.github.io/statmodels7/reference/path_steps.md)
  : The Step a Coordinate Descent Would Take on a Block
- [`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
  : The Values a Path Visits
- [`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
  : Factorize a Penalized Information Once
- [`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
  : The Log-Determinant of a Penalized Information, Robustly
- [`pd_logdet_dense()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet_dense.md)
  : The Dense Route of pd_logdet
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
- [`print.start_strategy`](https://statmodels7.github.io/statmodels7/reference/print.start_strategy.md)
  : Print a Starting-Value Strategy
- [`print_block()`](https://statmodels7.github.io/statmodels7/reference/print_block.md)
  : Print One Block of a Summary
- [`readable_hyper_rows()`](https://statmodels7.github.io/statmodels7/reference/readable_hyper_rows.md)
  : The Quantities a Penalty's Hyperparameters Are About
- [`refreshes_own_block()`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md)
  : Does a Term Recompute Its Own Block?
- [`reject_incompatible()`](https://statmodels7.github.io/statmodels7/reference/reject_incompatible.md)
  : Combinations of Terms That Are Not a Model
- [`reject_nested_offsets()`](https://statmodels7.github.io/statmodels7/reference/reject_nested_offsets.md)
  : Reject an Offset Buried Inside a Term
- [`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md)
  : Reject a Term the Fitting Scheme Does Not Cover
- [`row_nonzeros()`](https://statmodels7.github.io/statmodels7/reference/row_nonzeros.md)
  : A Design's Nonzeros, Ordered by Row
- [`search_coords()`](https://statmodels7.github.io/statmodels7/reference/search_coords.md)
  : Which Coefficients a Search Should Cover
- [`seg_grid_start()`](https://statmodels7.github.io/statmodels7/reference/seg_grid_start.md)
  : Choose a Break-Point Term's Starting Positions on a Grid
- [`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
  : The Smallest Admissible Value of a Shape Parameter
- [`short_keys()`](https://statmodels7.github.io/statmodels7/reference/short_keys.md)
  : A Term Key Shortened for Display
- [`smooth_linear_cols()`](https://statmodels7.github.io/statmodels7/reference/smooth_linear_cols.md)
  : Which Coefficients of a Smooth Are the Linear Part
- [`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
  : Invert a Matrix That Ought to Be Positive Definite
- [`sparse_augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/sparse_augmented_solve.md)
  : Solve a Scoring Step From a Sparse Square-Root Design
- [`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
  : Fit One Non-Smooth Block, the Others Held Fixed
- [`sparse_lmin()`](https://statmodels7.github.io/statmodels7/reference/sparse_lmin.md)
  : The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated
- [`spec_at()`](https://statmodels7.github.io/statmodels7/reference/spec_at.md)
  : Rebuild a Specification Against New Data
- [`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
  : Take the Offsets Out of an Equation
- [`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
  : The Square-Root Design
- [`start_at.StartIntercepts`](https://statmodels7.github.io/statmodels7/reference/start_at.StartIntercepts.md)
  : Starting Values From the Intercept-Only Fit
- [`start_at.StartOrigin`](https://statmodels7.github.io/statmodels7/reference/start_at.StartOrigin.md)
  : Starting Values at Zero
- [`start_at.StartRandom`](https://statmodels7.github.io/statmodels7/reference/start_at.StartRandom.md)
  : Starting Values From a Random Draw
- [`start_at.StartSearch`](https://statmodels7.github.io/statmodels7/reference/start_at.StartSearch.md)
  : Starting Values From a Global Search
- [`start_strategy_class()`](https://statmodels7.github.io/statmodels7/reference/start_strategy_class.md)
  : The start_strategy Class Object
- [`statmod_active()`](https://statmodels7.github.io/statmodels7/reference/statmod_active.md)
  : Which Coefficients Are Not Sitting at a Kink
- [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  : The Alternation Between the Smooth Block and the Rest
- [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  : Split a Specification's Terms Into the Smooth Block and the Rest
- [`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
  : Advance the Refresh State
- [`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
  : The Design at Given Coefficients
- [`statmod_design_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_blocks.md)
  : The Blocks of a Design
- [`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
  : Effective Degrees of Freedom, Per Term
- [`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
  : The Smoothing-Parameter Correction to the Effective Degrees of
  Freedom
- [`statmod_entry_key()`](https://statmodels7.github.io/statmodels7/reference/statmod_entry_key.md)
  : The Key of One of a Term's Penalties
- [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  : The Linear Predictors and the Parameters They Give
- [`statmod_filter_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
  : Run the Structural Terms at the Current Parameters
- [`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
  : Fit the Coefficients and a Filter's Parameters in One System
- [`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md)
  : Fit the Structural Terms' Own Parameters
- [`statmod_fitted_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_fitted_spec.md)
  : The Terms as the Fit Left Them
- [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
  : The Observed Information Over the Coefficients and a Filter's
  Parameters
- [`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md)
  : How Many Values a Path Visits for One Hyperparameter
- [`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
  : Which Hyperparameters the Terms Hold
- [`statmod_held_levels()`](https://statmodels7.github.io/statmodels7/reference/statmod_held_levels.md)
  : Which Structural Levels a Linear Intercept Already Carries
- [`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md)
  : Override the Starting Hyperparameters
- [`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md)
  : The Hyperparameters a Specification Starts From
- [`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)
  : The Variance of the Hyperparameters a Marginal Criterion Estimated
- [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
  : The Information of the Weighted Log-Likelihood
- [`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
  : The Intercept of Each Equation, on the Link Scale
- [`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
  : The Weighted Log-Likelihood of a Specification at Given Coefficients
- [`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
  : The Marginal Criterion at Given Coefficients and Hyperparameters
- [`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
  : The Penalized Curvature Over the Coefficients and a Filter's
  Parameters
- [`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
  : The Exact Gradient of the Marginal Criterion
- [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
  : The Exact Hessian of the Marginal Criterion
- [`statmod_min_ratio()`](https://statmodels7.github.io/statmodels7/reference/statmod_min_ratio.md)
  : How Far Down the Path Reaches for One Term
- [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  : The Objective, Its Gradient and Its Hessian, Stacked
- [`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
  : Select the Hyperparameters of a Kinked Penalty Along a Path
- [`statmod_path_setting()`](https://statmodels7.github.io/statmodels7/reference/statmod_path_setting.md)
  : One Setting of the Path, Read From the Term
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
- [`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md)
  : Have the Refreshable Terms Settled?
- [`statmod_refreshable()`](https://statmodels7.github.io/statmodels7/reference/statmod_refreshable.md)
  : Which Terms Recompute Their Own Block
- [`statmod_regime_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_at.md)
  : Run the Structural Terms at the Current Parameters
- [`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md)
  : The Observed Information of a Model Carrying a Regime Term
- [`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
  : The Same Model Read on Other Rows
- [`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
  : The Score of the Weighted Log-Likelihood
- [`statmod_search()`](https://statmodels7.github.io/statmodels7/reference/statmod_search.md)
  : How a Term Covers Its Own Hyperparameters
- [`statmod_select()`](https://statmodels7.github.io/statmodels7/reference/statmod_select.md)
  : Estimate the Hyperparameters, by Whichever Route Each One Admits
- [`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
  : Starting Coefficients
- [`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
  : Which Terms Rewrite the Likelihood
- [`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md)
  : The Exact Gradient Where a Penalty Covers a Filter's Own Parameters
- [`statmod_structural_par()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_par.md)
  : The Structural Terms' Estimated Parameters
- [`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
  : The Penalty Over a Structural Term's Own Parameters
- [`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)
  : The Score in a Structural Term's Own Parameters
- [`statmod_structural_state()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_state.md)
  : The State of the Structural Terms
- [`statmod_structural_table()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_table.md)
  : A Structural Term's Parameters, With Standard Errors
- [`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
  : Interpret and Build Each Parameter's Terms
- [`statmod_theta_shifted()`](https://statmodels7.github.io/statmodels7/reference/statmod_theta_shifted.md)
  : The Parameters Under One Regime
- [`statmod_unit()`](https://statmodels7.github.io/statmodels7/reference/statmod_unit.md)
  : One Penalized Unit, by Parameter and Key
- [`statmod_values()`](https://statmodels7.github.io/statmodels7/reference/statmod_values.md)
  : The Values a Term Wrote Out for One Hyperparameter
- [`statmodels7`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  [`statmodels7-package`](https://statmodels7.github.io/statmodels7/reference/statmodels7-package.md)
  : statmodels7: The S7 Toolkit for Statistical Modeling
- [`statmodels7_attach()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach.md)
  : Attach the Member Packages
- [`statmodels7_attach_message()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach_message.md)
  : The Body of the Attach Message
- [`structural_callbacks()`](https://statmodels7.github.io/statmodels7/reference/structural_callbacks.md)
  : The Score and Curvature a Filter Is Driven By
- [`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
  : The Two Pieces of the Chain Term That Read the Direction
- [`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md)
  : What the Joint Chain Term Needs Before a Direction Is Known
- [`structural_joint_basis()`](https://statmodels7.github.io/statmodels7/reference/structural_joint_basis.md)
  : The Subspace a Marginal Criterion Integrates Over, Jointly
- [`structural_kind()`](https://statmodels7.github.io/statmodels7/reference/structural_kind.md)
  : Which Shape of the Structural Contract a Term Implements
- [`structural_penalized()`](https://statmodels7.github.io/statmodels7/reference/structural_penalized.md)
  : Does a Structural Term Carry a Penalty of Its Own?
- [`structural_penalty_block()`](https://statmodels7.github.io/statmodels7/reference/structural_penalty_block.md)
  : The Penalty Over a Structural Term's Free Parameters, as a Block
- [`structural_psi()`](https://statmodels7.github.io/statmodels7/reference/structural_psi.md)
  : From the Unconstrained Scale to the Term's Parameters
- [`structural_range_cols()`](https://statmodels7.github.io/statmodels7/reference/structural_range_cols.md)
  : Which of a Structural Term's Free Parameters a Penalty Covers
- [`structural_zeta_start()`](https://statmodels7.github.io/statmodels7/reference/structural_zeta_start.md)
  : The Parameters a Structural Term Starts From
- [`summary(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  : Summarize a Fitted Model
- [`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
  : The Blocks of One Distribution Parameter
- [`term_block_kind()`](https://statmodels7.github.io/statmodels7/reference/term_block_kind.md)
  : What Kind of Block a Term Reports As
- [`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
  : Evaluate a Formula's Terms With modelterms7 in Front
- [`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
  : The Trace Against a Contraction, Without Forming It
- [`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
  : The Trace of the Determinant's Movement With the Coefficients
- [`unfittable_reason()`](https://statmodels7.github.io/statmodels7/reference/unfittable_reason.md)
  : Why a Term Is Outside the Fitting Scheme
- [`unknown_what()`](https://statmodels7.github.io/statmodels7/reference/unknown_what.md)
  : The Message for an Unrecognized Prediction Target
- [`vb_inner()`](https://statmodels7.github.io/statmodels7/reference/vb_inner.md)
  : The Verbosity of an Inner Fit Inside the Outer Search
- [`vb_name()`](https://statmodels7.github.io/statmodels7/reference/vb_name.md)
  : The Name of Whatever Is About to Run
- [`vb_rule()`](https://statmodels7.github.io/statmodels7/reference/vb_rule.md)
  : A Titled Rule for a Verbose Trace
- [`vb_say()`](https://statmodels7.github.io/statmodels7/reference/vb_say.md)
  : A Detail Line of a Verbose Trace
- [`vcov(`*`<StatmodFit>`*`)`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  : The Variance Matrix of a Fit
- [`verbosity()`](https://statmodels7.github.io/statmodels7/reference/verbosity.md)
  : Resolve the Verbosity Setting
- [`worth_sparse()`](https://statmodels7.github.io/statmodels7/reference/worth_sparse.md)
  : Is a Matrix Worth Factorizing Sparsely?
- [`zap_nonfinite()`](https://statmodels7.github.io/statmodels7/reference/zap_nonfinite.md)
  : Zero the Non-Finite Entries of a Penalty's Hessian
