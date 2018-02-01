__precompile__(true)

module QuadDIRECT

using StaticArrays, PositiveFactorizations
using Compat

export Box, leaves, splitprint, splitprint_red, treeprint, value, analyze, analyze!, minimize

include("types.jl")
include("util.jl")
include("algorithm.jl")

end # module
