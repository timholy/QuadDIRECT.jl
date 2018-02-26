function init(f, xsplits, lower::AbstractVector, upper::AbstractVector)
    # Validate the inputs
    n = length(lower)
    length(upper) == length(xsplits) == n || throw(DimensionMismatch("inconsistent dimensionality among lower, upper, and splits"))
    for i = 1:n
        xsi = xsplits[i]
        length(xsi) == 3 || throw(DimensionMismatch("3 values along each dimension required, got $(xsi) along dimension $i"))
        lower[i] <= xsi[1] < xsi[2] < xsi[3] <= upper[i] || error("splits must be in strictly increasing order and lie between lower and upper, got $xsi along dimension $i")
        all(isfinite, xsi) || error("not all splits along dimension $i are finite")
    end
    # Compute a consistent eltype
    xs1 = xsplits[1]
    Tx = promote_type(typeof(xs1[1]), typeof(xs1[2]), typeof(xs1[3]))
    for i = 2:n
        xsi = xsplits[i]
        Tx = promote_type(Tx, typeof(xsi[1]), typeof(xsi[2]), typeof(xsi[3]))
    end
    Tx = boxeltype(Tx)
    x0 = Tx[x[2] for x in xsplits]
    box, xstar = _init(f, copy(x0), xsplits, lower, upper)
    box, x0, xstar
end

@noinline function _init(f, xstar, xsplits, lower, upper)
    T = boxeltype(promote_type(eltype(xstar), eltype(lower), eltype(upper)))
    n = length(xstar)
    root = box = Box{T,n}()
    xtmp = copy(xstar)
    fcur = dummyvalue(T)
    for i = 1:n
        box = split!(box, f, xtmp, i, xsplits[i], lower[i], upper[i], xstar[i], fcur)
        xcur, fcur = box.parent.xvalues[box.parent_cindex], box.parent.fvalues[box.parent_cindex]
        xstar = replacecoordinate!(xstar, i, xcur)
    end
    box, xstar
end

function split!(box::Box{T}, f, xtmp, splitdim, xsplit, lower::Real, upper::Real, xcur, fcur) where T
    # Evaluate f along the splits, keeping track of the best
    fsplit = MVector3{T}(typemax(T), typemax(T), typemax(T))
    fmin, idxmin = typemax(T), 0
    for l = 1:3
        @assert(isfinite(xsplit[l]))
        if xsplit[l] == xcur && !isdummy(fcur)
            ftmp = fcur
        else
            xtmp = replacecoordinate!(xtmp, splitdim, xsplit[l])
            ftmp = f(xtmp)
        end
        if ftmp < fmin
            fmin = ftmp
            idxmin = l
        end
        fsplit[l] = ftmp
    end
    idxmin == 0 && error("function was not finite at any evaluation point $xsplit")
    if idxmin == 1 && fsplit[1] == fsplit[2]
        idxmin = 2  # prefer the middle in case of ties
    end
    add_children!(box, splitdim, xsplit, fsplit, lower, upper)
    xtmp = replacecoordinate!(xtmp, splitdim, xsplit[idxmin])
    return box.children[idxmin]
end

# Splits a box. If the "tilt" over the domain of the
# box suggests that one of its neighbors might be even lower,
# recursively calls itself to split that box, too.
# Returns `box, minval`.
function autosplit!(box::Box{T}, mes::Vector{<:MELink}, f::WrappedFunction, x0, xtmp, splitdim, xsplitdefaults, lower, upper, minwidth, visited::Set) where T
    box ∈ visited && error("already visited box")
    box.qnconverged && return (box, value(box))
    if !isleaf(box)
        # This box already got split along a different dimension
        box = find_smallest_child_leaf(box)
        xtmp = position(box, x0)
    end
    if splitdim == 0
        # If we entered this box from a neighbor, just choose an axis that has been split the least
        nsplits = count_splits(box)
        splitdim = indmin(nsplits)
        if nsplits[splitdim] > 0
            # If all dimensions have been split, prefer splitting unbounded dimensions
            bbs = boxbounds(box, lower, upper)
            if all(isfinite, bbs[splitdim])
                for (i, bb) in enumerate(bbs)
                    if !all(isfinite, bb)
                        splitdim = i
                        break
                    end
                end
            end
        end
    end
    xsplitdefault = xsplitdefaults[splitdim]
    lwr, upr = lower[splitdim], upper[splitdim]
    p = find_parent_with_splitdim(box, splitdim)
    if isroot(p) && p.splitdim != splitdim
        xcur, fcur = x0[splitdim], box.parent.fvalues[box.parent_cindex]
        split!(box, f, xtmp, splitdim, xsplitdefault, lwr, upr, xcur, fcur)
        trimschedule!(mes, box, splitdim, x0, lower, upper)
        return box, minimum(box.fvalues)
    end
    bb = boxbounds(p)
    bb[2]-bb[1] > max(minwidth[splitdim], epswidth(bb)) || return (box, value(box))
    xp, fp = p.parent.xvalues, p.parent.fvalues
    Δx = max(xp[2]-xp[1], xp[3]-xp[2])  # a measure of the "pragmatic" box scale even when the box is infinite
    xcur = xp[p.parent_cindex]
    @assert(xcur == xtmp[splitdim])
    fcur = box.parent.fvalues[box.parent_cindex]
    # Do a quadratic fit along splitdim
    xvert, fvert, qcoef = qfit(xp[1]=>fp[1], xp[2]=>fp[2], xp[3]=>fp[3])
    if qcoef > 0
        # xvert is a minimum. However, in general we can't trust it: for example, for the
        # "canyon" function
        #     f(x,y) = (x+y)^2/k + k*(x-y)^2,
        # at fixed `x` the minimum `y` is at
        #     ystar = (k^2-1)/(k^2+1)*x.
        # If the current box has a different `x` than the one used for box `p`,
        # then this estimate will be wrong.
        # (If `p == box.parent` it would be safe, but given that we'd have to pick a third
        # point anyway, this case does not seem worth treating separately.)
        # Instead of trusting `xvert`, we take `qcoef` (the quadratic coefficient) at face value. We then
        # choose a second (arbitrary) point along splitdim (e.g., a new `y` in the "canyon" example)
        # to evaluate `f`, and then use these two position/value pairs to estimate the
        # remaining two parameters of the quadratic.
        if isinf(bb[1])
            xnew = xcur - 8*(bb[2]-xcur)
        elseif isinf(bb[2])
            xnew = xcur + 8*(xcur-bb[1])
        else
            xnew = bb[2]-xcur > xcur-bb[1] ? (xcur+2*bb[2])/3 : (2*bb[1]+xcur)/3
        end
        @assert(isfinite(xnew))
        xtmp = replacecoordinate!(xtmp, splitdim, xnew)
        fnew = f(xtmp)
        # Solve for the minimum of the quadratic given the two pairs xcur=>fcur, xnew=>fnew
        xvert = (xcur+xnew)/2 - (fnew-fcur)/(2*qcoef*(xnew-xcur))
        if bb[1] <= xvert <= bb[2]
            # xvert is in the box, use it---but first, make sure it's not "degenerate" with
            # xcur or xnew
            xvert = ensure_distinct(xvert, xcur, xnew, bb)
            xtmp = replacecoordinate!(xtmp, splitdim, xvert)
            fvert = f(xtmp)
            xf1, xf2, xf3 = order_pairs(xcur=>fcur, xnew=>fnew, xvert=>fvert)
            add_children!(box, splitdim, MVector3{T}(xf1[1], xf2[1], xf3[1]),
                          MVector3{T}(xf1[2], xf2[2], xf3[2]), lwr, upr)
            trimschedule!(mes, box, splitdim, x0, lower, upper)
            return box, minimum(box.fvalues)
        end
        # xvert is not in the box. Prepare to split the neighbor, but for this box
        # just bisect xcur and xnew
        if xvert < bb[1]
            xtmp, dir = replacecoordinate!(xtmp, splitdim, bb[1]), -1
        else
            xtmp, dir = replacecoordinate!(xtmp, splitdim, bb[2]), +1
        end
        nbr, success = find_leaf_at_edge(get_root(box), xtmp, splitdim, dir)
        xmid = (xcur+xnew)/2
        xtmp = replacecoordinate!(xtmp, splitdim, xmid)
        fmid = f(xtmp)
        xf1, xf2, xf3 = order_pairs(xcur=>fcur, xnew=>fnew, xmid=>fmid)
        add_children!(box, splitdim, MVector3{T}(xf1[1], xf2[1], xf3[1]),
                      MVector3{T}(xf1[2], xf2[2], xf3[2]), lwr, upr)
        trimschedule!(mes, box, splitdim, x0, lower, upper)
        minbox = minimum(box.fvalues)
        if !success || nbr ∈ visited || !isfinite(value(nbr)) || length(visited) > ndims(box)
            return (box, minbox) # check nbr ∈ visited to avoid cycles
        end
        nbr, mv = autosplit!(nbr, mes, f, x0, position(nbr, x0), 0, xsplitdefaults, lower, upper, minwidth, push!(visited, box))
        if mv < minbox
            return nbr, mv
        end
        return box, minbox
    end
    trisect!(box, f, xtmp, splitdim, bb, Δx, xcur, fcur)
    trimschedule!(mes, box, splitdim, x0, lower, upper)
    return box, isleaf(box) ? value(box) : minimum(box.fvalues)
end

function trisect!(box::Box{T}, f, xtmp, splitdim, bb, Δx, xcur, fcur) where T
    l, r = bb
    if isinf(l)
        l = r - 6*Δx  # gap is usually 1/3 of the box size, so undo this at the level of "box size"
    elseif isinf(r)
        r = l + 6*Δx
    end
    a, b, c = 5l/6+r/6, (l+r)/2, l/6+5r/6
    # Ensure xcur is one of the three
    if xcur < (a+b)/2
        a = xcur
    elseif xcur < (b+c)/2
        b = xcur
    else
        c = xcur
    end
    if a < b < c
        return split!(box, f, xtmp, splitdim, MVector3{T}(a, b, c), bb..., xcur, fcur)
    end
    return box
end

# This requires there to be a (non-root) parent split along splitdim
function trisect!(box::Box, f, xtmp, splitdim, p = find_parent_with_splitdim(box, splitdim))
    bb = boxbounds(p)
    xp, fp = p.parent.xvalues, p.parent.fvalues
    Δx = max(xp[2]-xp[1], xp[3]-xp[2])  # a measure of the "pragmatic" box scale even when the box is infinite
    xcur = xp[p.parent_cindex]
    @assert(xcur == xtmp[splitdim])
    fcur = box.parent.fvalues[box.parent_cindex]
    trisect!(box, f, xtmp, splitdim, bb, Δx, xcur, fcur)
end

# Returns target box, value, and α. Uses α=0 to signal failure.
function quasinewton!(box::Box{T}, B, g, c, scale, f::Function, x0, lower, upper, itermax = 20) where T
    x = position(box, x0)
    fbox = fmin = value(box)
    isfinite(fbox) || return box, fmin, T(0)
    # Solve the bounded linear least-squares problem forcing B to be positive-definite
    cB = cholfact(Positive, B)
    Δx = -(cB \ g)   # unconstrained solution
    lsc, usc = (lower.-x)./scale, (upper.-x)./scale
    Δx[:] .= clamp.(Δx, lsc, usc)
    Δx = lls_bounded!(Δx, full(cB), g, lsc, usc)
    # Test whether the update moves the point appreciably. If not, mark this point
    # as converged.
    issame = true
    for i = 1:length(Δx)
        issame &= abs(Δx[i]) < sqrt(eps(scale[i]))
    end
    if issame
        box.qnconverged = true
        return box, fmin, T(0)
    end
    Δx[:] .= Δx .* scale

    α = T(1.0)
    iter = 0
    while !isinside(x + α*Δx, lower, upper) && iter < itermax
        α /= 2
        iter += 1
    end
    iter == itermax && return box, fmin, T(0)
    iter = 0 # the above weren't "real" iterations, so reset
    root = get_root(box)
    # Do a backtracking linesearch
    local xtarget, ftarget
    while iter < itermax
        iter += 1
        xtarget = x + α*Δx
        # Sadly, the following function evaluations get discarded. But recording them
        # while preserving the box structure would take ndims(box)-1 additional
        # evaluations, which is not worth it unless xtarget itself is an improvement.
        ftarget = f(xtarget)
        ftarget < fbox && break
        α /= 2
    end
    ftarget >= fbox && return box, fmin, T(0)

    leaf = find_leaf_at(root, xtarget)
    isfinite(value(leaf)) || return box, fmin, T(0)
    # The new point should improve on what was already obtained in `leaf`
    ftarget >= value(leaf) && return box, fmin, T(0)

    # For consistency with the tree structure, we can change only one coordinate of
    # xleaf per split. Cycle through all coordinates, picking the one at each stage
    # that yields the smallest value as predicted by the quadratic model among the
    # not-yet-selected coordinates.
    dims_targeted = falses(ndims(leaf))
    q(x) = (x'*B*x)/2 + g'*x + c
    for j = 1:ndims(leaf)
        xleaf = position(leaf, x0)
        xtest = copy(xleaf)
        imin, qmin = 0, typemax(T)
        for i = 1:ndims(leaf)
            dims_targeted[i] && continue
            xtest = replacecoordinate!(xtest, i, xtarget[i])
            qx = q(xtest)
            if qx < qmin
                qmin = qx
                imin = i
            end
        end
        imin == 0 && break
        xcur, xt = xleaf[imin], xtarget[imin]
        if abs(xcur - xt) < sqrt(eps(scale[imin]))
            # We're already at this point, continue
            dims_targeted[imin] = true
            continue
        end
        bb = boxbounds(leaf, imin, lower, upper)
        a, b, c = pick3(xcur, xt, bb)
        minsep = eps(T)*(c-a)
        if a + minsep >= b || b + minsep >= c
            # The box might be so narrow that there are not enough distinct floating-point
            # numbers that lie between the bounds.
            dims_targeted[imin] = true
            continue
        end
        split!(leaf, f, xleaf, imin, MVector3{T}(a, b, c), bb..., xleaf[imin], value(leaf))
        fmin = min(fmin, minimum(leaf.fvalues))
        childindex = findfirst(equalto(xt), leaf.xvalues)
        leaf = leaf.children[childindex]
        isfinite(value(leaf)) || return leaf, fmin, T(0)
        dims_targeted[imin] = true
    end
    return leaf, fmin, α
end

# Linear least squares with box-bounds. From
#    Sequential Coordinate-wise Algorithm for the
#    Non-negative Least Squares Problem
#    Vojtˇch e Franc, V ́clav aHlav ́ acˇ, and Mirko Navara
# Solves min(x'*H*x/2 + f'*x) over the domain bounded by lower and upper
function lls_bounded(H::AbstractMatrix{T}, f::VT, lower::V, upper::V, tol = sqrt(eps(T))) where {T,VT<:AbstractVector{T},V<:AbstractVector}
    x = clamp.(-H \ f, lower, upper)
    lls_bounded!(x, H, f, lower, upper, zeros(T, length(f)), tol)
    return x
end

lls_bounded!(x::VT, H::AbstractMatrix{T}, f::VT, lower::V, upper::V, tol = sqrt(eps(T))) where {T,VT<:AbstractVector{T},V<:AbstractVector} =
    lls_bounded!(x, H, f, lower, upper, zeros(T, length(f)), tol)

function lls_bounded!(x::VT, H::AbstractMatrix{T}, f::VT, lower::V, upper::V, g::VT, tol = sqrt(eps(T))) where {T,VT<:AbstractVector{T},V<:AbstractVector}
    assert_oneindexed(A) = @assert(all(r->first(r)==1, Compat.axes(A)))
    assert_oneindexed(x); assert_oneindexed(H); assert_oneindexed(f); assert_oneindexed(lower)
    assert_oneindexed(upper); assert_oneindexed(g)
    m, n = size(H)
    m == n || error("H must be square")
    length(f) == n || error("Mismatch between H and f")
    length(g) == n || error("Mismatch between H and g")
    niter = 0
    while niter <= 1000
        # Once per iter, calculate gradient freshly to avoid accumulation of roundoff
        A_mul_B!(g, H, x)
        for i = 1:n
            g[i] += f[i]
        end
        xchange = zero(T)
        for k = 1:n
            xold = x[k]
            xnew = clamp(xold - g[k]/H[k,k], lower[k], upper[k])
            x[k] = xnew
            dx = xnew-xold
            xchange = max(xchange, abs(dx))
            if dx != 0
                for i = 1:n
                    g[i] += dx*H[i,k]
                end
            end
        end
        if xchange < tol
            break
        end
        niter += 1
    end
    x
end

# A dumb O(N) algorithm for building the minimum-edge structures
# For large trees, this is the bottleneck
function minimum_edges(root::Box{T,N}, x0, lower, upper, minwidth=zeros(eltype(x0), ndims(root)); extrapolate::Bool=true) where {T,N}
    mes = [MELink{T,T}(root) for i = 1:N]
    # Make it a priority to have enough splits along each coordinate to provide adequate
    # information for the quasi-Newton method. To compute an estimate of the full quadratic
    # model, we'll need (N+1)*(N+2)÷2 independent points, so prioritze splitting
    # boxes so that each coordinate has at least 1/Nth of the requisite number of splits.
    nthresh = ceil(Int, (((N+1)*(N+2))÷2)/N)
    nsplits = Vector{Int}(uninitialized, N)
    for box in leaves(root)
        box.qnconverged && continue
        fval = value(box)
        isfinite(fval) || continue
        count_splits!(nsplits, box)
        nmin = minimum(nsplits)
        for i = 1:N
            if nmin < nthresh
                nsplits[i] == nmin || continue
            end
            bb = boxbounds(find_parent_with_splitdim(box, i), lower[i], upper[i])
            bb[2]-bb[1] < minwidth[i] && continue
            insert!(mes[i], width(box, i, x0, lower, upper), box=>extrapolate ? fval+qdelta(box, i) : fval)
        end
    end
    return mes
end

function trimschedule!(mes::Vector{<:MELink}, box::Box, splitdim, x0, lower, upper)
    isleaf(box) && return mes
    for child in box.children
        fval = child.parent.fvalues[child.parent_cindex]
        for i = 1:ndims(box)
            trim!(mes[i], width(child, i, x0, lower, upper), child=>fval)
        end
    end
    return mes
end

function mesprint(mes)
    for i = 1:length(mes)
        print(i, ": ")
        for item in mes[i]
            print(", (", item.w, ',', item.f, ')')
        end
        println()
    end
    nothing
end

function sweep!(root::Box, f, x0, splits, lower, upper; extrapolate::Bool = true, fvalue=-Inf, minwidth=zeros(eltype(x0), ndims(root)))
    mes = minimum_edges(root, x0, lower, upper, minwidth; extrapolate=extrapolate)
    sweep!(root, mes, f, x0, splits, lower, upper; fvalue=fvalue, minwidth=minwidth)
end
function sweep!(root::Box{T}, mes::Vector{<:MELink}, f, x0, splits, lower, upper; fvalue=-Inf, minwidth=zeros(eltype(x0), ndims(root))) where T
    xtmp = similar(x0)
    flag = similar(x0, Bool)
    nsplits = similar(x0, Int)
    visited = Set{typeof(root)}()
    splitboxes = Set{typeof(root)}()
    dimorder = sortperm(length.(mes))  # process the dimensions with the shortest queues first
    fvalueT = T(fvalue)
    for i in dimorder
        me = mes[i]
        while !isempty(me)
            item = popfirst!(me)
            box = item.l
            box.qnconverged && continue
            bb = boxbounds(find_parent_with_splitdim(box, i), lower[i], upper[i])
            bb[2]-bb[1] < minwidth[i] && continue
            position!(xtmp, flag, box)
            default_position!(xtmp, flag, x0)
            count_splits!(nsplits, box)
            if nsplits[i] > 0 && any(iszero, nsplits)
                # println("discarding ", box, " along dimension ", i)
                continue
            end
            empty!(visited)
            finalbox, minval = autosplit!(box, mes, f, x0, xtmp, i, splits, lower, upper, minwidth, visited)
            push!(splitboxes, finalbox)
            minval <= fvalueT && return root, unique_smallest_leaves(splitboxes)
        end
    end
    root, unique_smallest_leaves(splitboxes)
end

"""
    root, x0 = analyze(f, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)

Analyze the behavior of `f`, searching for minima, over the rectangular box specified by
`lower` and `upper` (`lower[i] <= x[i] <= upper[i]`). The bounds may be infinite.
`splits` is a list of 3-vectors containing the initial values along each coordinate
axis at which to potentially evaluate `f`; the values must be in increasing order.

`root` is a tree structure that stores information about the behavior of `f`.
`x0` contains the initial evaluation point, the position constructed from the middle value
of `splits` in each dimension.

`rtol` and `atol` represent relative and
absolute, respectively, changes in minimum function value required for the exploration
to terminate. (These limits must be hit on `ndims` successive sweeps.) Alternatively,
the analysis is terminated if the function value is ever reduced below `fvalue`.

# Example:

    # A function that has a long but skinny valley aligned at 45 degrees (minimum at [0,0])
    function canyon(x)
        x1, x2 = x[1], x[2]
        return 0.1*(x1+x2)^2 + 10*(x1-x2)^2
    end

    lower, upper = [-Inf, -Inf], [Inf, Inf]  # unbounded domain

    # We'll initially consider a grid that covers -11->-9 along the first dimension
    # and -7->-5 along the second. This isn't a great initial guess, but that's OK.
    splits = ([-11,-10,-9], [-7,-6,-5])

    # Since this problem has a minimum value of 0, relying on `rtol` is not ideal
    # (it takes quite a few iterations), so specify an absolute tolerance.
    julia> root, x0 = analyze(canyon, splits, lower, upper; atol=0.01)  # low-precision solution
    (BoxRoot@[NaN, NaN], [-10.0, -6.0])

    julia> box = minimum(root)
    Box0.015220406463743074@[0.192158, 0.19604]

    julia> value(box)
    0.015220406463743074

    julia> position(box, x0)
    2-element Array{Float64,1}:
    0.192158
    0.19604

    # Now a higher-precision solution
    julia> root, x0 = analyze(canyon, splits, lower, upper; atol=1e-5)
    (BoxRoot@[NaN, NaN], [-10.0, -6.0])

    julia> box = minimum(root)
    Box2.7058500379897107e-6@[-0.0025621, -0.00261386]

See also [`minimize`](@ref).
"""
function analyze(f, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500, print_interval=typemax(Int), kwargs...)
    box, x0 = init(f, splits, lower, upper)
    root = get_root(box)
    analyze!(root, f, x0, splits, lower, upper; rtol=rtol, atol=atol, fvalue=fvalue, maxevals=maxevals, print_interval=print_interval, kwargs...)
    return root, x0
end

"""
    root = analyze!(root, f, x0, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)

Further refinement of `root`. See [`analyze`](@ref) for details.
"""
function analyze!(root::Box, f::Function, x0, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500, print_interval=typemax(Int), kwargs...)
    analyze!(root, CountedFunction(f), x0, splits, lower, upper; rtol=rtol, atol=atol, fvalue=fvalue, maxevals=maxevals, print_interval=print_interval, kwargs...)
end
function analyze!(root::Box{T}, f::WrappedFunction, x0, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500, print_interval=typemax(Int), nquasinewton=qnthresh(ndims(root)), kwargs...) where T
    box = minimum(root)
    boxval = value(box)
    lastval = typemax(boxval)
    tol_counter = 0
    lastprint = 0
    baseline_evals = len = lenold = length(leaves(root))
    if baseline_evals == numevals(f)
        baseline_evals = 0   # already counted
    end
    if print_interval < typemax(Int)
        println("Initial minimum ($len evaluations): ", box)
    end
    extrapolate = true
    # The quasi-Newton step can reduce the function value so significantly, insist on
    # using it at least once (unless the threshold to kick in not reachable)
    used_quasinewton = nquasinewton > maxevals
    nsweep, nqn = baseline_evals+numevals(f), 0
    thresh = sqrt(eps(T))
    nsplits = zeros(Int, ndims(root))
    while boxval > fvalue && (tol_counter <= ndims(root) || !used_quasinewton) && len < maxevals
        lastval = boxval
        numeval0 = numevals(f)
        _, splitboxes = sweep!(root, f, x0, splits, lower, upper; extrapolate=extrapolate, fvalue=fvalue, kwargs...)
        box = minimum(root)
        numeval1 = numevals(f)
        nsweep += numeval1 - numeval0
        if nquasinewton < baseline_evals + numevals(f) < maxevals && boxval > fvalue && !isempty(splitboxes)
            # Quasi-Newton step for boxes split in the sweep
            spbxs = sort(collect(splitboxes); by=value) # if unsorted, set hashing would introduce randomness
            # determine which are in separate "basins" using a convexity test
            # basinboxes = spbxs
            # basinboxes = different_basins(spbxs, x0, lower, upper)
            basinboxes = filter(is_diag_convex, spbxs)
            badbox = Dict(box=>false for box in basinboxes)
            for (ibox, box) in enumerate(basinboxes)
                box.qnconverged && continue
                isleaf(box) || continue   # "outdated"
                badbox[box] && continue   # in same basin as a previous successful quadratic fit
                baseline_evals + numevals(f) >= maxevals && break
                # Ensure at least one split per coordinate before trying to build the quadratic model
                count_splits!(nsplits, box)
                xtmp = position(box, x0)
                for i = 1:ndims(box)
                    if nsplits[i] == 0
                        box = split!(box, f, xtmp, i, splits[i], lower[i], upper[i], xtmp[i], value(box))
                        xtmp[i] = box.parent.xvalues[box.parent_cindex]
                    end
                end
                scale = boxscale(box, splits)
                Q, xbase, c, success_build = build_quadratic_model(box, x0, scale, thresh)
                if success_build && all(x->x>0, Q.dimpiv)
                    # Add more points, as needed, to fill in any missing rows from the regression
                    @assert isleaf(box)
                    complete_quadratic_model!(Q, c, box, f, x0, xbase, scale, splits, lower, upper, thresh)
                    # If the model is complete, try Quasi-Newton algorithm
                    if Q.nzrows[] == size(Q.coefs, 1) && minimum(abs.(diag(Q.coefs))) >= thresh
                        g, B = solve(Q)
                        # Bscale = (1./scale) .* B .* (1./scale)'
                        # display(Bscale)
                        # error("stop")
                        leaf, minvalue, α = quasinewton!(box, B, g, c, scale, f, x0, lower, upper)
                        success_qn = α > 0
                        used_quasinewton |= success_qn
                        minvalue < fvalue && break
                        if success_qn
                            @assert(isleaf(leaf))
                            boxvalue = value(box)
                            # If the improvement was tiny, mark the solution to block
                            # future attempts at quasi-Newton in this region
                            if boxvalue - minvalue <= max(atol, rtol*(abs(boxvalue) + abs(minvalue)))
                                # Check that leaf has roughly the same coordinates, too.
                                xbox = position(box, x0)
                                xleaf = position(leaf, x0)
                                if issame(xbox, xleaf, scale)
                                    leaf.qnconverged = true
                                end
                            end
                            # check future boxes and eliminate them if they appear to lie in the same basin
                            if α == 1
                                for j = ibox+1:length(basinboxes)
                                    boxj = basinboxes[j]
                                    vj, xj = value(boxj), position(boxj, x0)
                                    Δx = (xj - xbase) ./ scale
                                    pred = c + g'*Δx + (Δx'*B*Δx)/2
                                    if abs(pred - vj) <= 0.1*abs(pred - minvalue)   # FIXME hardcoded
                                        badbox[boxj] = true
                                    end
                                end
                            end
                        end
                    end
                end
            end
            # if !used_quasinewton
            #     nquasinewton *= 2
            # end
            box = minimum(root)
            boxval = value(box)
        end
        nqn += numevals(f) - numeval1
        len = baseline_evals + numevals(f)
        len == lenold && break  # couldn't split any boxes
        lenold = len
        if len-lastprint > print_interval
            println("minimum ($len evaluations): ", box)
            lastprint = len
        end
        @assert(boxval <= lastval)
        if lastval - boxval < atol || lastval - boxval < rtol*(abs(lastval) + abs(boxval))
            tol_counter += 1
        else
            tol_counter = 0
        end
    end
    if print_interval < typemax(Int)
        println("Final minimum ($len evaluations): ", minimum(root))
        println("$nsweep evaluations for initialization+sweeps, $nqn for Quasi-Newton")
    end
    root
end

"""
    xmin, fmin = minimize(f, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)

Return the position `xmin` and value `fmin` of the minimum of `f` over the specified domain.
See [`analyze`](@ref) for information about the input arguments.
"""
function minimize(f, splits, lower, upper; kwargs...)
    root, x0 = analyze(f, splits, lower, upper; kwargs...)
    box = minimum(root)
    return position(box, x0), value(box)
end
