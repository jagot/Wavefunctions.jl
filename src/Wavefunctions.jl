module Wavefunctions

using LinearMaps
using RecipesBase

import Base: norm, normalize!, *, A_mul_B!

abstract type AbstractWfn{T} end

# Single-active electron wavefunctions

# An AbstractSAEWfn is assumed to have a data matrix, containing the
# wavefunction with increasing spatial dimension along the rows, and
# increasing partial wave index along the columns.
abstract type AbstractSAEWfn{T} <: AbstractWfn{T} end

function normalize!(wfn::AbstractSAEWfn)
    scale!(wfn.data, 1/norm(wfn))
    wfn
end

function A_mul_B!(Y::W, A::Union{AbstractMatrix,AbstractLinearMap}, B::W) where W <: AbstractSAEWfn
    A_mul_B!(view(Y.data, :), A, view(B.data, :))
    Y
end

include("wfn_utils.jl")
include("eigenstates.jl")
include("fd_sae_wfn.jl")

export AbstractWfn

end # module
