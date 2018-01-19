"""
   xvert, fvert, qcoef = qfit(xm=>fm, x0=>f0, xp=>fp)

Given three points `xm < x0 < xp ` and three corresponding
values `fm`, `f0`, and `fp`, fit a quadratic. Returns the position `xvert` of the vertex,
the quadratic's value `fvert` at `xvert`, and the coefficient `qcoef` of the quadratic term.
`xvert` is a minimum if `qcoef > 0`.

Note that if the three points lie on a line, `qcoef == 0` and both `xvert` and `fvert` will
be infinite.
"""
function qfit(xfm, xf0, xfp)
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    @assert(xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp))
    cm = fm/((xm-x0)*(xm-xp))  # coefficients of Lagrange polynomial
    c0 = f0/((x0-xm)*(x0-xp))
    cp = fp/((xp-xm)*(xp-x0))
    qcoef = cm+c0+cp
    qvalue(x) = cm*(x-x0)*(x-xp) + c0*(x-xm)*(x-xp) + cp*(x-xm)*(x-x0)
    if fm == f0 == fp
        return x0, f0, zero(qcoef)
    end
    xvert = (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*qcoef)
    return xvert, qvalue(xvert), qcoef
end

## Minimum Edge List utilities
Base.empty!(mel::MELink) = (mel.next = mel; return mel)

function dropnext!(prev, next)
    if next != next.next
        # Drop the next item from the list
        next = next.next
        prev.next = next
    else
        # Drop the last item from the list
        prev.next = prev
        next = prev
    end
    return next
end

function Base.insert!(mel::MELink, w, lf::Pair)
    l, f = lf
    prev, next = mel, mel.next
    while prev != next && w > next.w
        if f <= next.f
            next = dropnext!(prev, next)
        else
            prev = next
            next = next.next
        end
    end
    if w == next.w
        f >= next.f && return mel
        next = dropnext!(prev, next)
    end
    if prev == next
        # we're at the end of the list
        prev.next = typeof(mel)(l, w, f)
    else
        if f < next.f
            prev.next = typeof(mel)(l, w, f, next)
        end
    end
    return mel
end

Base.start(mel::MELink) = mel
Base.done(mel::MELink, state::MELink) = state == state.next
Base.next(mel::MELink, state::MELink) = (state.next, state.next)

function Base.show(io::IO, mel::MELink)
    print(io, "List(")
    next = mel.next
    while mel != next
        print(io, '(', next.w, ", ", next.l, "=>", next.f, "), ")
        mel = next
        next = next.next
    end
    print(io, ')')
end


## Box utilities
function Base.show(io::IO, box::Box)
    x = fill(NaN, ndims(box))
    position!(x, box)
    val = isroot(box) ? "Root" : value(box)
    print(io, "Box$val@", x)
end

function value(box::Box)
    isroot(box) && error("root box does not have a unique value")
    box.parent.fvalues[box.parent_cindex]
end
value_safe(box::Box{T}) where T = isroot(box) ? typemax(T) : value(box)

Base.isless(box1::Box, box2::Box) = isless(value_safe(box1), value_safe(box2))

function treeprint(io::IO, f::Function, root::Box)
    show(io, root)
    y = f(root)
    y != nothing && print(io, y)
    if !isleaf(root)
        print(io, '(')
        treeprint(io, f, root.children[1])
        print(io, ", ")
        treeprint(io, f, root.children[2])
        print(io, ", ")
        treeprint(io, f, root.children[3])
        print(io, ')')
    end
end
treeprint(io::IO, root::Box) = treeprint(io, x->nothing, root)

function add_children!(parent::Box, splitdim, xvalues, fvalues, u::Real, v::Real)
    isleaf(parent) || error("cannot add children to non-leaf node")
    (length(xvalues) == 3 && xvalues[1] < xvalues[2] < xvalues[3]) || throw(ArgumentError("xvalues must be monotonic, got $xvalues"))
    parent.splitdim = splitdim
    p = find_parent_with_splitdim(parent, splitdim)
    if isroot(p)
        parent.minmax = (u, v)
    else
        parent.minmax = boxbounds(p)
    end
    parent.xvalues = xvalues
    parent.fvalues = fvalues
    for i = 1:3
        Box(parent, i)  # creates the children of parent
    end
    parent
end

function cycle_free(box)
    p = parent(box)
    while !isroot(p)
        p == box && return false
        p = p.parent
    end
    return true
end

function isparent(parent, child)
    parent == child && return true
    while !isroot(child)
        child = child.parent
        parent == child && return true
    end
    return false
end

"""
    boxp = find_parent_with_splitdim(box, splitdim::Integer)

Return the first node at or above `box` who's parent box was split
along dimension `splitdim`.
"""
function find_parent_with_splitdim(box::Box, splitdim::Integer)
    while !isroot(box)
        p = parent(box)
        if p.splitdim == splitdim
            return box
        end
        box = p
    end
    return box
end

"""
    box = find_smallest_child_leaf(root)

Walk the tree recursively, choosing the child with smallest function value at each stage.
`box` will be a leaf node.
"""
function find_smallest_child_leaf(box::Box)
    # Not guaranteed to be the smallest function value, it's the smallest that can be
    # reached stepwise
    while !isleaf(box)
        idx = indmin(box.fvalues)
        box = box.children[idx]
    end
    box
end

"""
    box = find_leaf_at_edge(root, x, splitdim, dir)

Return the leaf-node `box` that contains `x` with an edge at `x[splitdim]`.
If `dir > 0`, a box to the right of the edge will be returned; if `dir < 0`, a box to the
left will be returned.

This is a useful utility for finding the neighbor of a given box. Example:

    # Let's find the neighbors of `box` along its parent's splitdim
    x = position(box, x0)
    i = box.parent.splitdim
    bb = boxbounds(box)
    # Right neighbor
    x[i] = bb[2]
    rnbr = find_leaf_at_edge(root, x, i, +1)
    # Left neighbor
    x[i] = bb[1]
    lnbr = find_leaf_at_edge(root, x, i, -1)
"""
function find_leaf_at_edge(root::Box, x, splitdim::Integer, dir::Signed)
    isleaf(root) && return root
    while !isleaf(root)
        i = root.splitdim
        found = false
        for box in (root.children[1], root.children[2], root.children[3])
            bb = boxbounds(box)
            if within(x[i], bb, dir)
                root = box
                found = true
                break
            end
        end
        found || error("$(x[i]) not within $(root.minmax)")
    end
    root
end

"""
    x = position(box)
    x = position(box, x0)

Return the n-dimensional position vector `x` at which this box was evaluated
when it was a leaf. Some entries of `x` might be `NaN`, if `box` is sufficiently
near the root and not all dimensions have been split.
The variant supplying `x0` fills in those dimensions with the corresponding values
from `x0`.
"""
position(box::Box) = position!(fill(NaN, ndims(box)), box)

function position(box::Box, x0::AbstractVector)
    x = fill(NaN, ndims(box))
    flag = falses(length(x0))
    position!(x, flag, box)
    default_position!(x, flag, x0)
end

function position!(x, box::Box)
    flag = falses(length(x))
    position!(x, flag, box)
    return x
end
function position!(x, flag, box::Box)
    fill!(flag, false)
    nfilled = 0
    while !isroot(box) && nfilled < length(x)
        i = box.parent.splitdim
        if !flag[i]
            x[i] = box.parent.xvalues[box.parent_cindex]
            flag[i] = true
            nfilled += 1
        end
        box = box.parent
    end
    x
end
function default_position!(x, flag, xdefault)
    length(x) == length(flag) == length(xdefault) || throw(DimensionMismatch("all three inputs must have the same length"))
    for i = 1:length(x)
        if !flag[i]
            x[i] = xdefault[i]
        end
    end
    x
end

"""
    left, right = boxbounds(box)

Compute the bounds of `box` along the `splitdim` of `box`'s parent.
This throws an error for the root box.
"""
function boxbounds(box::Box)
    isroot(box) && error("cannot compute bounds on root Box")
    p = parent(box)
    if box.parent_cindex == 1
        return (p.minmax[1], (p.xvalues[1]+p.xvalues[2])/2)
    elseif box.parent_cindex == 2
        return ((p.xvalues[1]+p.xvalues[2])/2, (p.xvalues[2]+p.xvalues[3])/2)
    elseif box.parent_cindex == 3
        return ((p.xvalues[2]+p.xvalues[3])/2, p.minmax[2])
    end
    error("invalid parent_cindex $(box.parent_cindex)")
end

"""
    left, right = boxbounds(box, lower::Real, upper::Real)

Compute the bounds of `box` along the `splitdim` of `box`'s parent.
For the root box, returns `(lower, upper)`.
"""
function boxbounds(box::Box, lower::Real, upper::Real)
    isroot(box) && return (lower, upper)
    return boxbounds(box)
end

"""
    bb = boxbounds(box, lower::AbstractVector, upper::AbstractVector)

Compute the bounds of `box` along all dimensions.
"""
function boxbounds(box::Box{T}, lower::AbstractVector, upper::AbstractVector) where T
    length(lower) == length(upper) == ndims(box) || throw(DimensionMismatch("lower and upper must match dimensions of box"))
    bb = [(T(lower[i]), T(upper[i])) for i = 1:ndims(box)]
    boxbounds!(bb, box)
end

function boxbounds!(bb, box::Box)
    flag = falses(ndims(box))
    boxbounds!(bb, flag, box)
    return bb
end
function boxbounds!(bb, flag, box::Box)
    fill!(flag, false)
    if isleaf(box)
        bb[box.parent.splitdim] = boxbounds(box)
        flag[box.parent.splitdim] = true
    else
        bb[box.splitdim] = box.minmax
        flag[box.splitdim] = true
    end
    nfilled = 1
    while !isroot(box) && nfilled < ndims(box)
        i = box.parent.splitdim
        if !flag[i]
            bb[i] = boxbounds(box)
            flag[i] = true
            nfilled += 1
        end
        box = box.parent
    end
    bb
end

function width(box::Box, splitdim::Integer, xdefault::Real, lower::Real, upper::Real)
    p = find_parent_with_splitdim(box, splitdim)
    bb = boxbounds(p, lower, upper)
    x = isroot(p) ? xdefault : p.parent.xvalues[p.parent_cindex]
    max(x-bb[1], bb[2]-x)
end
width(box::Box, splitdim::Integer, xdefault, lower, upper) =
    width(box, splitdim, xdefault[splitdim], lower[splitdim], upper[splitdim])

"""
    within(x, (left, right), dir)

Return `true` if `x` lies between `left` and `right`. If `x` is on the edge,
`dir` must point towards the interior (positive if `x==left`, negative if `x==right`).
"""
function within(x::Real, bb::Tuple{Real,Real}, dir)
    if !(bb[1] <= x <= bb[2])
        return false
    end
    ((x == bb[1]) & (dir < 0)) && return false
    ((x == bb[2]) & (dir > 0)) && return false
    return true
end

function Base.extrema(root::Box)
    isleaf(root) && error("tree is empty")
    minv, maxv = extrema(root.fvalues)
    for bx in root
        isleaf(bx) && continue
        mn, mx = extrema(bx.fvalues)
        minv = min(minv, mn)
        maxv = max(maxv, mx)
    end
    minv, maxv
end

## Tree traversal
function get_root(box::Box)
    while !isroot(box)
        box = parent(box)
    end
    box
end

abstract type DepthFirstIterator end
Base.iteratorsize(::Type{<:DepthFirstIterator}) = Base.SizeUnknown()

struct DepthFirstLeafIterator{T} <: DepthFirstIterator
    root::Box{T}
end

struct VisitorBool{T}
    box::Box{T}
    done::Bool
end

function visit_leaves(root::Box)
    DepthFirstLeafIterator(root)
end

function Base.start(iter::DepthFirstLeafIterator)
    find_next_leaf(iter, VisitorBool(iter.root, false))
end
Base.start(root::Box) = VisitorBool(root, false)

Base.done(iter::DepthFirstLeafIterator, state::VisitorBool) = state.done
Base.done(root::Box, state::VisitorBool) = state.done

function Base.next(iter::DepthFirstLeafIterator, state::VisitorBool)
    @assert(isleaf(state.box))
    return (state.box, find_next_leaf(iter, state))
end
function find_next_leaf(iter::DepthFirstLeafIterator, state::VisitorBool)
    _, state = next(iter.root, state)
    while !isleaf(state.box) && !state.done
        _, state = next(iter.root, state)
    end
    return state
end

function Base.next(root::Box, state::VisitorBool)
    item, done = state.box, state.done
    if isleaf(item)
        box, i = up(item, root)
        if i <= length(box.children)
            return (item, VisitorBool(box.children[i], false))
        end
        @assert(box == root)
        return (item, VisitorBool(root, true))
    end
    return (item, VisitorBool(item.children[1], false))
end

function up(box, root)
    local i
    while true
        box, i = box.parent, box.parent_cindex+1
        box == root && return (box, i)
        i <= length(box.children) && break
    end
    return (box, i)
end

## Utilities for working with both mutable and immutable vectors
replacecoordinate!(x, i::Integer, val) = (x[i] = val; x)

replacecoordinate!(x::SVector{N,T}, i::Integer, val) where {N,T} =
    SVector{N,T}(_rpc(Tuple(x), i-1, T(val)))
@inline _rpc(t, i, val) = (ifelse(i == 0, val, t[1]), _rpc(tail(t), i-1, val)...)
_rps(::Tuple{}, i, val) = ()

ipcopy!(dest, src) = copy!(dest, src)
ipcopy!(dest::SVector, src) = src

## Other utilities
lohi(x, y) = x <= y ? (x, y) : (y, x)
function lohi(x, y, z)
    @assert(x <= y)
    z <= x && return z, x, y
    z <= y && return x, z, y
    return x, y, z
end

function biggest_interval(a, b, c, d)
    ab, bc, cd = b-a, c-b, d-c
    if ab <= bc && ab <= cd
        return (a, b)
    elseif bc <= ab && bc <= cd
        return (b, c)
    end
    return (c, d)
end
