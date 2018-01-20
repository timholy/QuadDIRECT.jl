module QuadDIRECT

using StaticArrays
using Compat

export leaves, treeprint, value, analyze, minimize

include("types.jl")
include("util.jl")
include("algorithm.jl")

end # module
