module ReversePropagation

export gradient, forward_backward_contractor, SSAFunction

import Symbolics: toexpr, variable

using SymbolicUtils
using SymbolicUtils: Sym, Term
using SymbolicUtils.Rewriters

using Symbolics
using Symbolics: value,
                iscall, operation, arguments,
                Assignment

using IntervalContractors
using IntervalArithmetic, IntervalArithmetic.Symbols

import IntervalArithmetic.Symbols: ⊓, ⊔


@register_symbolic a ⊓ b false
@register_symbolic a ⊔ b false

using OrderedCollections

using ChainRulesCore, ChainRules

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)



# struct Assignment
#     lhs
#     rhs
# end

# Base.show(io::IO, eq::Assignment) = print(io, lhs(eq), " := ", rhs(eq))

include("make_variable.jl")
include("ssa.jl")
include("chain_rules_interface.jl")
include("cse.jl")
include("reverse_diff.jl")
include("reverse_icp.jl")

end
