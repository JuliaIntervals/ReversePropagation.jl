# ReversePropagation.jl release notes

## v0.6.0

### Breaking

- `ChainRules` and `ChainRulesCore` are no longer dependencies. Scalar
  derivative rules live in an in-package dispatch-based table
  (`src/scalar_rules.jl`) with `adj(f, z̄, xs...)` as the public
  reverse-mode interface. The external behaviour of `gradient` /
  `forward_backward_contractor` is unchanged; downstream code that
  called `ReversePropagation.adj` via a `ChainRules.rrule`-style
  signature will need to adapt.
- `max` / `min` moved from the unary to the binary op registry (they
  are genuinely binary — the previous unary listing was a bug). They
  are now fully supported end-to-end, with subgradient-aware rules.
- `sign` is no longer in the reversibility registry either. Its
  derivative is distributional at 0 and zero elsewhere, so it has no
  well-defined place in a pipeline whose whole point is symbolic
  differentiation. Expressions containing `sign(...)` now fail at
  gradient/contractor build time with a clear error rather than
  silently producing a mathematically meaningless result.
  `_abs_subgrad`'s `Real`-input branch still calls `Base.sign` at
  runtime internally, but it never appears in any SSA.

### Added

- Subgradient-aware scalar rules for `abs`, `max`, `min`. Evaluate to
  ordinary derivatives for scalar `Real` inputs; for `Interval` /
  `BareInterval` inputs at (or across) a kink, return the subgradient
  interval — e.g. `∂|x|/∂x` over `[-2, 3]` → `[-1, 1]`.
- Type-preserving `log(2)` / `log(10)` helpers (`_ln2`, `_ln10`) so the
  adjoints of `exp2` / `exp10` / `log2` / `log10` stay rigorous when
  evaluated over intervals. Previously the constant was snapped to a
  `Float64` literal.
- Exported `binarize_ssa(ssa::SSAFunction)` — decomposes n-ary or
  compound-argument assignments (e.g. `_ā := _b̄ * cos(_a)`) into
  elementary unary/binary ops. Useful as a post-processing step on
  the result of `gradient_code` before feeding into downstream
  interval-reverse consumers.
- Regression test that sweeps every differentiable entry in the
  registries and asserts that `gradient_code |> binarize_ssa` produces
  SSA whose rhs's are all in the reversibility registry or the set of
  declared derivative-rule helpers (`DIFF_HELPERS`).

## v0.5.1

### Changed

- Runtime code generation in `gradient` and `forward_backward_contractor`
  now uses `@RuntimeGeneratedFunction` instead of `eval`. The user-facing
  calling convention is unchanged, but callers can now build and invoke a
  generated function within the same dynamic extent without hitting Julia
  world-age errors. Build-time is also significantly faster (4–15× on the
  benchmarked cases); steady-state call cost is unchanged.

### Added

- `RuntimeGeneratedFunctions` as a direct dependency.

## v0.5.0

### Breaking

- `cse_equations`, `forward_backward_code`, and `gradient_code` now return an
  `SSAFunction` instead of a `(code, final, …)` tuple. Callers that previously
  destructured the tuple should instead access `.code`, `.output`, and
  `.variables` (a `NamedTuple` carrying extras such as `constraint` or
  `gradient`).

### Added

- New exported `SSAFunction` type representing a symbolic function in SSA
  (single static assignment) form: a vector of `Assignment`s plus an output
  variable and a `NamedTuple` of additional named variables.
- `forward_backward_code`, `forward_backward_expr`, `forward_backward_contractor`,
  and `gradient_code` each gain an `SSAFunction`-accepting method, so a
  pre-computed SSA (e.g. the output of `gradient_code`) can be fed directly
  into `forward_backward_contractor`. The expression-accepting methods become
  thin wrappers that call `cse_equations` first.
- `show` method for `SSAFunction` that pretty-prints the assignments, and
  displays reverse operations as `rev(+)(…)` rather than raw generated names.

### Fixed

- Stale `istree` import (renamed to `iscall` in newer SymbolicUtils).
