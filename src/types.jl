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
function Base.convert(::Type{MVector3{T}}, a::AbstractVector) where T
    @boundscheck indices(a) == (Base.OneTo(3),) || throw(DimensionMismatch("vector must have length 3"))
    @inbounds ret = MVector3{T}(a[1], a[2], a[3])
    return ret
end
Base.convert(::Type{MVector3}, a::AbstractVector{T}) where T =
    convert(MVector3{T}, a)

Base.similar(::MVector3{T}, ::Type{S}, dims::Tuple{Vararg{Int}}) where {T,S} =
    Array{S}(uninitialized, dims)

boxtype(::Type{<:Integer}) = Float64
boxtype(::Type{T}) where T = T

# TODO
# 1. Define the Box structure
# 2. Define traversal routines for extracting 3-pt approximations along coordinate axes
# 3. Define split routines
# 4. Define outer API (iteration and termination criteria)
# 5. Consider collecting enough points to build the full quadratic model
mutable struct Box{T,N}
    parent::Box{T,N}
    parent_cindex::Int    # of its parent's children, which one is this?
    splitdim::Int         # the dimension along which this box has been split, or 0 if this is a leaf node
    minmax::Tuple{T,T}    # the outer edges not corresponding to one of the parent's xvalues (splits that occur between dots in Fig 2)
    xvalues::MVector3{T}  # the values of x_splitdim at which f is evaluated
    fvalues::MVector3{T}  # the corresponding values of f
    children::MVector3{Box{T,N}}

    function default!(box::Box{T,N}) where {T,N}
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

