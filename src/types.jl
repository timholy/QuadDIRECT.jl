"""
    mel = MELink{Tx,Ty}(dummylabel)

Create a linked-list storing labels and x/y pairs that represent the "minimum edge" of a
set of values, where the stored x and y values are monotonic in both x and y. Insertion
of a new x/y pair expunges any entries with smaller x but larger y. Add new elements with
`insert!(mel, x, label=>y)`.

The first item in the list is a dummy item with `x=0`, so that the head of the list is constant.
"""
mutable struct MELink{Tl,Tw,Tf}
    l::Tl   # label
    w::Tw   # x-coordinate (box width)
    f::Tf   # y-coordinate (function value)
    next::MELink{Tl,Tw,Tf}

    function MELink{Tl,Tw,Tf}(l, w, f) where {Tw, Tf, Tl}
        mel = new{Tl,Tw,Tf}(l, w, f)
        mel.next = mel
        return mel
    end
    function MELink{Tl,Tw,Tf}(l, w, f, next) where {Tw, Tf, Tl}
        mel = new{Tl,Tw,Tf}(l, w, f, next)
        return mel
    end
end

MELink{Tw,Tf}(dummylabel::Tl) where {Tl,Tw,Tf} = MELink{Tl,Tw,Tf}(dummylabel, typemin(Tw), typemin(Tf))

Base.iteratorsize(::Type{<:MELink}) = Base.SizeUnknown()

# Mutable length-3 vector without the storage overhead of Array
mutable struct MVector3{T} <: AbstractVector{T}
    x1::T
    x2::T
    x3::T
end
Base.IndexStyle(::Type{<:MVector3}) = IndexLinear()
Base.size(v::MVector3) = (3,)
Base.indices1(v::MVector3) = Base.OneTo(3)
function Base.getindex(v::MVector3, i::Int)
    @boundscheck ((1 <= i) & (i <= 3)) || Base.throw_boundserror(v, i)
    ifelse(i == 1, v.x1, ifelse(i == 2, v.x2, v.x3))
end
function Base.setindex!(v::MVector3, val, i::Int)
    @boundscheck ((1 <= i) & (i <= 3)) || Base.throw_boundserror(v, i)
    if i == 1
        v.x1 = val
    elseif i == 2
        v.x2 = val
    else
        v.x3 = val
    end
    v
end
Base.convert(::Type{MVector3{T}}, a::MVector3{T}) where T = a
function Base.convert(::Type{MVector3{T}}, a::AbstractVector) where T
    @boundscheck indices(a) == (Base.OneTo(3),) || throw(DimensionMismatch("vector must have length 3"))
    @inbounds ret = MVector3{T}(a[1], a[2], a[3])
    return ret
end
Base.convert(::Type{MVector3}, a::AbstractVector{T}) where T =
    convert(MVector3{T}, a)

Base.similar(::MVector3{T}, ::Type{S}, dims::Tuple{Vararg{Int}}) where {T,S} =
    Array{S}(uninitialized, dims)

boxeltype(::Type{<:Integer}) = Float64
boxeltype(::Type{T}) where T = T

mutable struct Box{T,N}
    parent::Box{T,N}
    parent_cindex::Int    # of its parent's children, which one is this?
    splitdim::Int         # the dimension along which this box has been split, or 0 if this is a leaf node
    qnconverged::Bool     # true if the box was the minimum-value box in a Quasi-Newton step and the improvement was within convergence criterion
    minmax::Tuple{T,T}    # the outer edges not corresponding to one of the parent's xvalues (splits that occur between dots in Fig 2)
    xvalues::MVector3{T}  # the values of x_splitdim at which f is evaluated
    fvalues::MVector3{T}  # the corresponding values of f
    children::MVector3{Box{T,N}}

    function default!(box::Box{T,N}) where {T,N}
        box.qnconverged = false
        box.minmax = (typemax(T), typemin(T))
        box.xvalues = MVector3{T}(typemax(T), zero(T), typemin(T))
        box.fvalues = MVector3{T}(zero(T), zero(T), zero(T))
        box.children = MVector3(box, box, box)
    end
    function Box{T,N}() where {T,N}
        # Create the root box
        box = new{T,N}()
        box.parent = box
        box.parent_cindex = 0
        box.splitdim = 0
        default!(box)
        return box
    end
    function Box{T,N}(parent::Box, parent_cindex::Integer) where {T,N}
        # Create a new child and store it in the parent
        box = new{T,N}(parent, parent_cindex, 0)
        parent.children[parent_cindex] = box
        default!(box)
        box
    end
end

Box(parent::Box{T,N}, parent_cindex::Integer) where {T,N} = Box{T,N}(parent, parent_cindex)

isroot(box::Box) = box.parent == box
isleaf(box::Box) = box.splitdim == 0
Base.parent(box::Box) = box.parent
Base.ndims(box::Box{T,N}) where {T,N} = N
Base.iteratorsize(::Type{<:Box}) = Base.SizeUnknown()

abstract type WrappedFunction <: Function end

# A function that keeps track of the number of evaluations
mutable struct CountedFunction{F} <: WrappedFunction
    f::F
    evals::Int

    CountedFunction{F}(f) where F = new{F}(f, 0)
end
CountedFunction(f::Function) = CountedFunction{typeof(f)}(f)

function (f::CountedFunction{F})(x::AbstractVector) where F
    f.evals += 1
    return f.f(x)
end

numevals(f::CountedFunction) = f.evals

# A function that keeps track of the number of evaluations
mutable struct LoggedFunction{F} <: WrappedFunction
    f::F
    values::Vector{Float64}

    LoggedFunction{F}(f) where F = new{F}(f, Float64[])
end
LoggedFunction(f::Function) = LoggedFunction{typeof(f)}(f)

function (f::LoggedFunction{F})(x::AbstractVector) where F
    val = f.f(x)
    push!(f.values, val)
    return val
end

numevals(f::LoggedFunction) = length(f.values)

# Quadratic-model Incremental Gaussian Elimination
struct QmIGE{T,N}
    coefs::PermutedDimsArray{T,2,(2, 1),(2, 1),Array{T,2}} # to make row-major operations fast
    rhs::Vector{T}
    dimpiv::Vector{Int}       # order in which dimensions were added
    rowzero::Vector{Bool}     # true if a row in coefs is all-zeros
    rowbox::Vector{Box{T,N}}  # the box that set this row
    ndims::typeof(Ref(1))     # number of dimensions added so far
    nzrows::typeof(Ref(1))
    rowtmp::Vector{T}
end

function QmIGE{T,N}() where {T,N}
    m = ((N+1)*(N+2))÷2 - 1  # the constant is handled separately
    coefs = PermutedDimsArray(zeros(T, m, m), (2, 1))
    rhs = zeros(T, m)
    rowtmp = zeros(T, m)
    dimpiv = zeros(Int, N)
    rowzero = fill(true, m)
    rowbox = Vector{Box{T,N}}(uninitialized, m)
    QmIGE{T,N}(coefs, rhs, dimpiv, rowzero, rowbox, Ref(0), Ref(0), rowtmp)
end
