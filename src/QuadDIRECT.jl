module QuadDIRECT

using StaticArrays, PositiveFactorizations
using Compat

export leaves, treeprint, value, analyze, analyze!, minimize

include("types.jl")
include("util.jl")
include("algorithm.jl")

end # module
