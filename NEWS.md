# ReversePropagation.jl release notes

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
