dummyvalue(::Type{T}) where T<:AbstractFloat = T(NaN)
dummyvalue(::Type{T}) where T = typemax(T)

isdummy(val::T) where T = isequal(val, dummyvalue(T))

"""
    qnthresh(N)

Return the minimum number of points needed to specify the quasi-Newton quadratic model
in `N` dimensions.
"""
qnthresh(N) = ((N+1)*(N+2))÷2

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
    cm, c0, cp = lagrangecoefs(xfm, xf0, xfp)
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    qvalue(x) = cm*(x-x0)*(x-xp) + c0*(x-xm)*(x-xp) + cp*(x-xm)*(x-x0)
    qcoef = cm+c0+cp
    if fm == f0 == fp
        return x0, f0, zero(qcoef) # when it's flat, use the middle point as the "vertex"
    end
    xvert = (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*qcoef)
    return xvert, qvalue(xvert), qcoef
end

@inline function lagrangecoefs(xfm, xf0, xfp)
    xm, fm = xfm.first, xfm.second
    x0, f0 = xf0.first, xf0.second
    xp, fp = xfp.first, xfp.second
    @assert(xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp))
    cm = fm/((xm-x0)*(xm-xp))  # coefficients of Lagrange polynomial
    c0 = f0/((x0-xm)*(x0-xp))
    cp = fp/((xp-xm)*(xp-x0))
    cm, c0, cp
end

"""
    Δf = qdelta(box)

Return the difference `fmin - fbox`, where `fbox` is the value of `f` at the evaluation
point of `box` and `fmin` is the minimum of a one-dimensional quadratic fit of `f` (using
the data in `box.parent`) over the bounds of `box`.

Note that `Δf` might be -Inf, if `box` is unbounded and the quadratic estimate
is not convex.
"""
function qdelta(box::Box)
    xv, fv = box.parent.xvalues, box.parent.fvalues
    xm, x0, xp = xv[1], xv[2], xv[3]
    fm, f0, fp = fv[1], fv[2], fv[3]
    fbox = box.parent.fvalues[box.parent_cindex]
    cm, c0, cp = lagrangecoefs(xm=>fm, x0=>f0, xp=>fp)
    qvalue(x) = cm*(x-x0)*(x-xp) + c0*(x-xm)*(x-xp) + cp*(x-xm)*(x-x0)
    bb = boxbounds(box)
    qcoef = cm+c0+cp
    xvert = (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*qcoef)
    if qcoef > 0 && bb[1] <= xvert <= bb[2]
        # Convex and the vertex is inside the box
        return qvalue(xvert) - fbox
    end
    # Otherwise, the minimum is achieved at one of the edges
    # This needs careful evaluation in the case of infinite boxes to avoid Inf - Inf == NaN.
    if isinf(bb[1]) || isinf(bb[2])
        lcoef = -cm*(x0+xp) - c0*(xm+xp) - cp*(xm+x0)
        ccoef = cm*x0*xp + c0*xm*xp + cp*xm*x0
        if qcoef == 0  # the function is linear
            let lcoef = lcoef, ccoef = ccoef
                lvalue(x) = lcoef*x + ccoef
                return min(lvalue(bb[1]), lvalue(bb[2])) - fbox
            end
        else
            qvalue_inf(x) = isinf(x) ? abs(x)*sign(qcoef) : qvalue(x)
            return min(qvalue_inf(bb[1]), qvalue_inf(bb[2])) - fbox
        end
    end
    return min(qvalue(bb[1]), qvalue(bb[2])) - fbox
end

function qdelta(box::Box{T}, splitdim::Integer) where T
    p = find_parent_with_splitdim(box, splitdim)
    return p.parent.splitdim == splitdim ? qdelta(p) : zero(T)
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

"""
    prev, next = trim!(mel::MELink, w, label=>fvalue)

Remove entries from `mel` that are worse than `(w, label=>fvalue)`.
Returns linked-list positions on either side of the putative new value.
"""
function trim!(mel::MELink, w, lf::Pair)
    l, f = lf.first, lf.second
    prev, next = mel, mel.next
    while prev != next && w > next.w
        if f <= next.f
            next = dropnext!(prev, next)
        else
            prev = next
            next = next.next
        end
    end
    if w == next.w && f < next.f
        next = dropnext!(prev, next)
    end
    return prev, next
end

function Base.insert!(mel::MELink, w, lf::Pair)
    l, f = lf
    prev, next = trim!(mel, w, lf)
    # Perform the insertion
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

function popfirst!(mel::MELink)
    item = mel.next
    mel.next = item.next == item ? mel : item.next
    return item
end

Base.length(mel::MELink) = count(x->true, mel)

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


### Box utilities
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

function pick_other(xvalues, fvalues, idx)
    j = 1
    if j == idx j += 1 end
    xf1 = xvalues[j] => fvalues[j]
    j += 1
    if j == idx j += 1 end
    xf2 = xvalues[j] => fvalues[j]
    return xf1, xf2
end

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
along dimension `splitdim`. If `box` has not yet been split along
`splitdim`, returns the root box.
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
    box = greedy_smallest_child_leaf(root)

Walk the tree recursively, choosing the child with smallest function value at each stage.
`box` will be a leaf node.
"""
function greedy_smallest_child_leaf(box::Box)
    # Not guaranteed to be the smallest function value, it's the smallest that can be
    # reached stepwise
    while !isleaf(box)
        idx = indmin(box.fvalues)
        box = box.children[idx]
    end
    box
end

function find_smallest_child_leaf(box::Box)
    vmin, boxmin = value(box), box
    for leaf in leaves(box)
        vleaf = value(leaf)
        if vleaf <= vmin
            vmin, boxmin = vleaf, leaf
        end
    end
    return boxmin
end

function unique_smallest_leaves(boxes)
    uboxes = Set{eltype(boxes)}()
    for box in boxes
        push!(uboxes, find_smallest_child_leaf(box))
    end
    return uboxes
end

"""
    box = find_leaf_at(root, x)

Return the leaf-node `box` that contains `x`.
"""
function find_leaf_at(root::Box, x)
    isleaf(root) && return root
    while !isleaf(root)
        i = root.splitdim
        found = false
        for box in (root.children[1], root.children[2], root.children[3])
            bb = boxbounds(box)
            if bb[1] <= x[i] <= bb[2]
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
    box, success = find_leaf_at_edge(root, x, splitdim, dir)

Return the node `box` that contains `x` with an edge at `x[splitdim]`.
If `dir > 0`, a box to the right of the edge will be returned; if `dir < 0`, a box to the
left will be returned. `box` will be a leaf-node unless `success` is false;
the usual explanation for this is that the chosen edge is at the edge of the allowed
domain, and hence there isn't a node in that direction.

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
    isleaf(root) && return root, false
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
        found || return (root, false)
    end
    return (root, true)
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
Base.position(box::Box) = position!(fill(NaN, ndims(box)), box)

function Base.position(box::Box, x0::AbstractVector)
    x = similar(x0)
    flag = falses(length(x0))
    position!(x, flag, box, x0)
end
function position!(x, flag, box::Box, x0::AbstractVector)
    copy!(x, x0)
    position!(x, flag, box)
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

"""
    bb = boxbounds(box, splitdim::Integer, lower::AbstractVector, upper::AbstractVector)

Compute the bounds of `box` along dimension `splitdim`.
"""
function boxbounds(box::Box{T}, splitdim::Integer, lower::AbstractVector, upper::AbstractVector) where T
    p = find_parent_with_splitdim(box, splitdim)
    boxbounds(p, lower[splitdim], upper[splitdim])
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
function boxbounds!(bb, flag, box::Box, lower, upper)
    for i = 1:ndims(box)
        bb[i] = (lower[i], upper[i])
    end
    boxbounds!(bb, flag, box)
end

function width(box::Box, splitdim::Integer, xdefault::Real, lower::Real, upper::Real)
    p = find_parent_with_splitdim(box, splitdim)
    bb = boxbounds(p, lower, upper)
    x = isroot(p) ? xdefault : p.parent.xvalues[p.parent_cindex]
    max(x-bb[1], bb[2]-x)
end
width(box::Box, splitdim::Integer, xdefault, lower, upper) =
    width(box, splitdim, xdefault[splitdim], lower[splitdim], upper[splitdim])

function isinside(x, lower, upper)
    ret = true
    for i = 1:length(x)
        ret &= lower[i] <= x[i] <= upper[i]
    end
    ret
end
function isinside(x, bb::Vector{Tuple{T,T}}) where T
    ret = true
    for i = 1:length(x)
        bbi = bb[i]
        ret &= bbi[1] <= x[i] <= bbi[2]
    end
    ret
end

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

function count_splits(box::Box)
    nsplits = Vector{Int}(uninitialized, ndims(box))
    count_splits!(nsplits, box)
end

function count_splits!(nsplits, box::Box)
    fill!(nsplits, 0)
    if box.splitdim == 0
        box = box.parent
        box.splitdim == 0 && return nsplits  # an unsplit tree
    end
    nsplits[box.splitdim] += 1
    while !isroot(box)
        box = box.parent
        nsplits[box.splitdim] += 1
    end
    return nsplits
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

## Utilities for experimenting with topology of the tree
function splitprint(io::IO, box::Box)
    if isleaf(box)
        print(io, 'l')
    else
        print(io, box.splitdim, '(')
        splitprint(io, box.children[1])
        print(io, ", ")
        splitprint(io, box.children[2])
        print(io, ", ")
        splitprint(io, box.children[3])
        print(io, ')')
    end
end
splitprint(box::Box) = splitprint(STDOUT, box)

function splitprint_red(io::IO, box::Box, thisbox::Box)
    if isleaf(box)
        box == thisbox ? print_with_color(:light_red, io, 'l') : print(io, 'l')
    else
        box == thisbox ? print_with_color(:light_red, io, box.splitdim) : print(io, box.splitdim)
        print(io, '(')
        splitprint_red(io, box.children[1], thisbox)
        print(io, ", ")
        splitprint_red(io, box.children[2], thisbox)
        print(io, ", ")
        splitprint_red(io, box.children[3], thisbox)
        print(io, ')')
    end
end
splitprint_red(box::Box, thisbox::Box) = splitprint_red(STDOUT, box, thisbox)

function Base.parse(::Type{B}, str::AbstractString) where B<:Box
    b = B()
    splitbox!(b, str)
end

# splitbox! uses integer-valued positions and sets all function values to 0
function splitbox!(box::Box{T,N}, dim) where {T,N}
    x = position(box, zeros(N))
    xd = x[dim]
    add_children!(box, dim, [xd,xd+1,xd+2], zeros(3), -Inf, Inf)
    box
end

function splitbox!(box::Box, str::AbstractString)
    str == "l" && return
    m = match(r"([0-9]*)\((.*)\)", str)
    dim = parse(Int, m.captures[1])
    dimstr = m.captures[2]
    splitbox!(box, dim)
    commapos = [0,0]
    commaidx = 0
    open = 0
    i = start(dimstr)
    while !done(dimstr, i)
        c, i = next(dimstr, i)
        if c == '('
            open += 1
        elseif c == ')'
            open -= 1
        elseif c == ',' && open == 0
            commapos[commaidx+=1] = prevind(dimstr, i)
        end
    end
    splitbox!(box.children[1], strip(dimstr[1:prevind(dimstr, commapos[1])]))
    splitbox!(box.children[2], strip(dimstr[nextind(dimstr, commapos[1]):prevind(dimstr, commapos[2])]))
    splitbox!(box.children[3], strip(dimstr[nextind(dimstr, commapos[2]):end]))
    return box
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

struct DepthFirstLeafIterator{B<:Box} <: DepthFirstIterator
    root::B
end

struct VisitorBool{B<:Box}
    box::B
    done::Bool
end

function leaves(root::Box)
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
    box == root && return (box, length(box.children)+1)
    while true
        box, i = box.parent, box.parent_cindex+1
        box == root && return (box, i)
        i <= length(box.children) && break
    end
    return (box, i)
end

function Base.length(iter::DepthFirstLeafIterator)
    state = start(iter)
    len = 0
    while !done(iter, state)
        _, state = next(iter, state)
        len += 1
    end
    len
end

## Utilities for working with both mutable and immutable vectors
replacecoordinate!(x, i::Integer, val) = (x[i] = val; x)

replacecoordinate!(x::SVector{N,T}, i::Integer, val) where {N,T} =
    SVector{N,T}(_rpc(Tuple(x), i-1, T(val)))
@inline _rpc(t, i, val) = (ifelse(i == 0, val, t[1]), _rpc(Base.tail(t), i-1, val)...)
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

function order_pairs(xf1, xf2)
    x1, f1 = xf1
    x2, f2 = xf2
    return x1 <= x2 ? (xf1, xf2) : (xf2, xf1)
end

function order_pairs(xf1, xf2, xf3)
    xf1, xf2 = order_pairs(xf1, xf2)
    xf2, xf3 = order_pairs(xf2, xf3)
    xf1, xf2 = order_pairs(xf1, xf2)
    return xf1, xf2, xf3
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

function ensure_distinct(x::T, xref, bb::Tuple{Real,Real}; minfrac = 0.1) where T
    Δx = min(xref - bb[1], bb[2] - xref)
    if !isfinite(Δx)
        x != xref && return x
        return xref + 1
    end
    Δxmin = T(minfrac*Δx)
    if abs(x - xref) < Δxmin
        s = x == xref ? (bb[2] - xref > xref - bb[1] ? 1 : -1) : sign(x-xref)
        x = T(xref + Δxmin*s)
    end
    return x
end

function ensure_distinct(x::T, x1, x2, bb::Tuple{Real,Real}; minfrac = 0.1) where T
    x1, x2 = lohi(x1, x2)
    Δxmin = minfrac*(x2-x1)
    if abs(x - x1) < Δxmin
        s = x == x1 ? 1 : sign(x-x1)
        x = T(max(bb[1], x1 + Δxmin*s))
    elseif abs(x - x2) < Δxmin
        s = x == x2 ? -1 : sign(x-x2)
        x = T(min(bb[2], x2 + Δxmin*s))
    end
    return x
end

"""
    t, exitdim = pathlength_box_exit(x0, dx, bb)

Given a ray `x0 + t*dx`, compute the value of `t` at which the ray exits the box
delimited by `bb` (a vector of `(lo, hi)` tuples). Also return the coordinate dimension
along which the exit occurs.
"""
function pathlength_box_exit(x0, dx, bb)
    t = oftype((bb[1][1] - x0[1])/dx[1], Inf)
    exitdim = 0
    for i = 1:length(x0)
        xi, dxi, bbi = x0[i], dx[i], bb[i]
        ti = (ifelse(dxi >= 0, bbi[2], bbi[1]) - xi)/dxi
        if ti < t
            t = ti
            exitdim = i
        end
    end
    t, exitdim
end

"""
    t, intersectdim = pathlength_hyperplane_intersect(x0, dx, xtarget, tmax)

Compute the maximum pathlength `t` (up to a value of `tmax`) at which the ray `x0 + t*dx`
intersects one of the hyperplanes specified by `x[i] = xtarget[i]` for any dimension `i`.
"""
function pathlength_hyperplane_intersect(x0, dx, xtarget, tmax)
    t = zero(typeof((xtarget[1] - x0[1])/dx[1]))
    intersectdim = 0
    for i = 1:length(x0)
        xi, dxi, xti = x0[i], dx[i], xtarget[i]
        ti = (xti - xi)/dxi
        if ti <= tmax && ti > t
            t = ti
            intersectdim = i
        end
    end
    t, intersectdim
end

function different_basins(boxes::AbstractVector{B}, x0, lower, upper) where B<:Box
    basinboxes = B[]
    for box in boxes
        ubasin = true
        for bbox in basinboxes
            if !is_different_basin(box, bbox, x0, lower, upper)
                ubasin = false
                break
            end
        end
        if ubasin
            push!(basinboxes, box)
        end
    end
    return basinboxes
end

function is_different_basin(box1, box2, x0, lower, upper)
    # This convexity test has its limits: we're comparing the maximum along the secant
    # within the box to a function value estimated from a diagonal quadratic model.
    # False positives happen.
    root = get_root(box1)
    v1, v2 = value(box1), value(box2)
    x1, x2 = position(box1, x0), position(box2, x0)
    dx = x2 - x1
    bb = boxbounds(box1, lower, upper)
    flag = Vector{Bool}(uninitialized, length(lower))
    leaf = box1
    t, exitdim = pathlength_box_exit(x1, dx, bb)
    # For consistency (e.g., commutivity with box1 and box2), we have to check
    # exit condition even at first box, even though it's likely a false positive.
    vmax = v1 + t*(v2-v1)
    qdtot = oftype(vmax, 0)
    for i = 1:ndims(leaf)
        qdtot += qdelta(leaf, i)
    end
    if value(leaf)+qdtot > vmax
        return true
    end
    while t < 1
        x2[:] .= x1 .+ t.*dx
        x2[exitdim] = bb[exitdim][dx[exitdim] > 0 ? 2 : 1]  # avoid roundoff error in the critical coordinate
        leaf_old = leaf
        leaf, success = find_leaf_at_edge(root, x2, exitdim, dx[exitdim] > 0 ? +1 : -1)
        success || break
        boxbounds!(bb, flag, leaf, lower, upper)
        tnext, exitdim = pathlength_box_exit(x1, dx, bb)
        tnext = max(tnext, t)
        while tnext == t  # must have hit a corner
            tnext += oftype(t, 1e-4)
            x2[:] .= x1 .+ tnext.*dx
            bxtmp = find_leaf_at(root, x2)
            boxbounds!(bb, flag, bxtmp, lower, upper)
            tnext, exitdim = pathlength_box_exit(x1, dx, bb)
        end
        vmax = v1 + t*(v2-v1)
        if tnext < 1
            vmax = max(vmax, v1 + tnext*(v2-v1))
        end
        qdtot = oftype(vmax, 0)
        for i = 1:ndims(leaf)
            qdtot += qdelta(leaf, i)
        end
        if value(leaf)+qdtot > vmax
            return true
        end
        t = tnext
    end
    return false
end

"""
    a, b, c = pick3(a, b, (lower::Real, upper::Real))

Returns an ordered triple `a, b, c`, with two agreeing with the input `a` and `b`,
and the third point bisecting the largest interval between `a`, `b`, and the edges
`lower`, `upper`.
"""
function pick3(a, b, bb)
    a, b = lohi(a, b)
    imin, imax = biggest_interval(bb[1], a, b, bb[2])
    if isinf(imin)
        return a-2*(b-a), a, b
    elseif isinf(imax)
        return a, b, b+2*(b-a)
    end
    return a, b, c = lohi(a, b, (imin+imax)/2)
end
