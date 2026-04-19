# ReversePropagation.jl release notes

## v0.6.0

### Breaking

- `ChainRules` and `ChainRulesCore` are no longer dependencies. Scalar
  derivative rules live in an in-package table (`src/scalar_rules.jl`)
  with `adj(f, z̄, xs...)` as the public reverse-mode interface. The
  external behaviour of `gradient` / `forward_backward_contractor` is
  unchanged; downstream code that called `ReversePropagation.adj` via
  a `ChainRules.rrule`-style signature will need to adapt.
- `sign`, `max`, `min` have been removed from the unary op registry.
  `sign`'s derivative is distributional and has no reverse contractor;
  `max` / `min` are binary and require interval-subgradient reverses
  that aren't in place yet. Expressions that used these in the forward
  pipeline will now error explicitly instead of silently building
  contractors that can't run.

### Added

- Exported `binarize_ssa(ssa::SSAFunction)` — decomposes n-ary or
  compound-argument assignments (e.g. `_ā := _b̄ * cos(_a)`) into
  elementary unary/binary ops. Useful as a post-processing step on
  the result of `gradient_code` before feeding into downstream
  interval-reverse consumers.
- Regression test that sweeps every entry in the unary/binary registries
  and asserts that `gradient_code |> binarize_ssa` produces SSA whose
  rhs's are all registered (reversible) operations.

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
