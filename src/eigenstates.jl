struct Eigenstate{T<:Real,W<:AbstractWfn{T}}
    wfn::W
    E::T
end

function Eigenstate(ee::LinAlg.Eigen, n::Int, WfnType::Type{<:AbstractWfn}, args...)
    n in eachindex(ee[:values]) || error("Requested eigenstate of invalid n quantum number")
    wfn = WfnType(vec2mat(Complex, ee[:vectors][:,n]), args...)
    normalize!(wfn)
    Eigenstate(wfn, ee[:values][n])
end

Eigenstate(H::AbstractMatrix, n, WfnType, args...) =
    Eigenstate(eigfact(H), n, WfnType, args...)

function Eigenstate(H::AbstractLinearMap, n::Int, WfnType::Type{<:AbstractWfn}, args...)
    E, φ = eigs(H, which = :SR, nev = n)
    wfn = WfnType(vec2mat(Complex, φ[:,n]), args...)
    normalize!(wfn)
    Eigenstate(wfn, E[n])
end

export Eigenstate
