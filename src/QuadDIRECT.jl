__precompile__(true)

module QuadDIRECT

using StaticArrays, PositiveFactorizations
using Compat

export Box, CountedFunction, LoggedFunction
export leaves, value, numevals, analyze, analyze!, minimize
export splitprint, splitprint_red, treeprint

include("types.jl")
include("util.jl")
include("algorithm.jl")

# Convenience methods
function analyze(f, x0::AbstractVector{<:Real}, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}; kwargs...)
    inds = linearindices(x0)
    linearindices(lower) == linearindices(upper) == inds || error("lengths must match")
    splits = similar(x0, Vector{Float64})
    for i in inds
        lo, hi, xi = lower[i], upper[i], x0[i]
        lo = max(xi-1, lo)
        hi = min(xi+1, hi)
        if lo == xi || hi == xi
            xi = (lo+hi)/2
        end
        @assert(lo < xi < hi)
        @assert isfinite(lo) && isfinite(xi) && isfinite(hi)
        splits[i] = [lo, xi, hi]
    end
    return analyze(f, splits, lower, upper; kwargs...)
end

function analyze(f, x0::AbstractVector{<:Real}; kwargs...)
    return analyze(f, x0, fill(-Inf, length(x0)), fill(Inf, length(x0)); kwargs...)
end

function analyze(f, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}; kwargs...)
    x0 = (lower .+ upper)./2
    for i in eachindex(x0)
        if !isfinite(x0[i])
            x0[i] = clamp(0, lower[i], upper[i])
        end
    end
    return analyze(f, x0, lower, upper; kwargs...)
end

end # module
