# Finite-difference, single-active electron
struct FDSAEWfn{T} <: AbstractSAEWfn{T}
    data::Matrix{Complex{T}}
    dr::T
end

norm(wfn::FDSAEWfn) = vecnorm(wfn.data)*sqrt(wfn.dr)

*(A::Union{AbstractMatrix,AbstractLinearMap}, wfn::FDSAEWfn) =
    FDSAEWfn(reshape(A*view(wfn.data, :), size(wfn.data)), wfn.dr)

@recipe function plot(wfn::FDSAEWfn)
    wfn.dr*(1:size(wfn.data,1)), abs2.(wfn.data)
end

export FDSAEWfn
