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
    Tx = boxtype(Tx)
    x0 = Tx[x[2] for x in xsplits]
    box, xstar = _init(f, copy(x0), xsplits, lower, upper)
    box, x0, xstar
end

@noinline function _init(f, xstar, xsplits, lower, upper)
    T = boxtype(promote_type(eltype(xstar), eltype(lower), eltype(upper)))
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
    idxmin == 0 && error("function was not finite at any evaluation point")
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
function autosplit!(box::Box{T}, mes::Vector{<:MELink}, f, x0, xtmp, splitdim, xsplitdefaults, lower, upper, minwidth, visited::Set) where T
    box ∈ visited && error("already visited box")
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
        return box
    end
    bb = boxbounds(p)
    bb[2]-bb[1] >= minwidth[splitdim] || return box
    N = ndims(box)
    nqn = ((N+1)*(N+2))÷2  # number of points needed for quasi-Newton approach
    if npoints_exceeds(box, nqn)
        coefs, values, positions, success = gather_independent_points(box, x0, nqn)
        if success
            success = quasinewton!(box, mes, coefs, values, f, x0, splitdim, lower, upper)
            success && return box
        end
    end
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
            return box
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
        (!success || nbr ∈ visited) && return box # don't get into a cycle
        return autosplit!(nbr, mes, f, x0, position(nbr, x0), 0, xsplitdefaults, lower, upper, minwidth, push!(visited, box))
    end
    # Trisect
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
    split!(box, f, xtmp, splitdim, MVector3{T}(a, b, c), bb..., xcur, fcur)
    trimschedule!(mes, box, splitdim, x0, lower, upper)
    return box
end

# Use regression to compute the best-fit quadratic model (with a dense Hessian)
function full_quadratic_fit(coefs::AbstractMatrix{T}, values::AbstractVector{T}) where T
    sz = size(coefs, 1)
    N = ceil(Int, sqrt(2*sz))
    while (N+1)*(N+2)÷2 > sz
        N -= 1
    end
    @assert((N+1)*(N+2) == 2*size(coefs, 1))
    B = Matrix{T}(uninitialized, N, N)
    g = Vector{T}(uninitialized, N)
    params = svdfact(coefs') \ values
    k = 0
    for j = 1:N
        B[j, j] = params[k+=1]
        for i = 1:j-1
            B[i, j] = B[j, i] = params[k+=1]
        end
    end
    for i = 1:N
        g[i] = params[k+=1]
    end
    c = params[end]
    return B, g, c
end

function setcol!(coefs, col, x, xref)
    N = length(x)
    nparams = ((N+1)*(N+2))÷2
    k = (col - 1)*nparams
    if length(coefs) < k+nparams
        resize!(coefs, k+nparams)
    end
    # Coefficients of the Hessian
    for j = 1:N
        dxj = x[j] - xref[j]
        coefs[k+=1] = 0.5*dxj*dxj
        for i = 1:j-1
            coefs[k+=1] = (x[i]-xref[i])*dxj
        end
    end
    # Coefficients of the gradient
    for i = 1:N
        coefs[k+=1] = x[i]-xref[i]
    end
    # Coefficient of the constant
    coefs[k+=1] = 1
    coefs
end

function quasinewton!(box::Box{T}, mes, coefs, values, f, x0, splitdim, lower, upper, itermax = 20) where T
    B, g, c = full_quadratic_fit(coefs, values)
    cB = cholfact(Positive, B)
    Δx = -(cB \ g)
    α = T(1.0)
    x = position(box, x0)
    iter = 0
    while !isinside(x + α*Δx, lower, upper) && iter < itermax
        α /= 2
        iter += 1
    end
    iter == itermax && return false
    iter = 0 # the above weren't "real" iterations, so reset
    root = get_root(box)
    # Do a backtracking linesearch, splitting any boxes that we encounter along the way
    while iter < itermax
        iter += 1
        xtarget = x + α*Δx
        # TODO: rewrite this to try f(xtarget), and if it's lower start splitting boxes.
        # (If not, backtrack without making any boxes. This will break the correspondence btw
        # function evaluation and leaves. Bummer.)
        # The problem with the current implementation is that even a great step triggers
        # backtracking because the coordinate-aligned approximation is too coarse.
        leaf = find_leaf_at(root, xtarget)
        xleaf = position(leaf, x0)
        # If leaf or one of its ancestors has been targeted before from an "external" box,
        # terminate. The only allowed re-targetings are from inside the narrowest box yet
        # targeted. This prevents running many line searches that all point to the same minimum.
        if leaf != box
            dims_targeted = qtargeted(leaf, x, lower, upper)
            if all(dims_targeted)
                iter == 1 && record_targeted!(mes, leaf, splitdim)
                return false
            end
        else
            dims_targeted = falses(ndims(leaf))
        end
        leaf.qtargeted = true
        # For consistency with the tree structure, we can change only one coordinate of
        # xleaf. Pick the one that yields the smallest value as predicted by the
        # quadratic model.
        xtest = copy(xleaf)
        q(x) = (x'*B*x)/2 + g'*x + c
        imin, fmin = 0, typemax(T)
        for i = 1:ndims(leaf)
            dims_targeted[i] && continue
            xtest = replacecoordinate!(xtest, i, xtarget[i])
            qx = q(xtest)
            if qx < fmin
                fmin = qx
                imin = i
            end
        end
        bb = boxbounds(leaf, imin, lower, upper)
        xcur = xleaf[imin]
        xt = ensure_distinct(xtarget[imin], xcur, bb)
        a, b, c = pick3(xcur, xt, bb)
        split!(leaf, f, xleaf, imin, MVector3{T}(a, b, c), bb..., xleaf[imin], value(leaf))
        # return true
        # end
        # # We're going to evaluate f inside leaf. For consistency with the tree structure,
        # # we need to match the position of the evaluation point in leaf at all but one
        # # coordinate. First, match the coordinate with the farthest intersection between
        # # the ray and the coordinate hyperplanes.
        # bbs = boxbounds(leaf, lower, upper)
        # texit = pathlength_box_exit(x, Δx, bbs)[1]
        # @show texit
        # tmax = min(oftype(texit, 1), texit)
        # t, intersectdim = pathlength_hyperplane_intersect(x, Δx, xleaf, tmax)
        # # At t, it matches xleaf[intersectdim]. We need to match N-2 remaining coordinates.
        # @assert(t > 0 && intersectdim != 0)
        # @show t intersectdim
        # xtarget = x + t*Δx
        # @show x xtarget xleaf
        # # Pick a 3rd point
        # bb = bbs[intersectdim]
        # a, b, c = pick3(xleaf[intersectdim], xtarget[intersectdim], bb)
        # split!(leaf, f, xleaf, intersectdim, MVector3{T}(a, b, c), bb..., xleaf[intersectdim], value(leaf))
        if minimum(leaf.fvalues) < value(box) || leaf == box
            return true
        end
        α /= 2
    end
    return false
end

# Might want to add leaf to mes?
record_targeted!(mes, leaf, splitdim) = nothing

# A dumb O(N) algorithm for building the minimum-edge structures
# For large trees, this is the bottleneck
function minimum_edges(root::Box{T,N}, x0, lower, upper, minwidth=zeros(eltype(x0), ndims(root)); extrapolate::Bool=true) where {T,N}
    mes = [MELink{T,T}(root) for i = 1:N]
    for box in leaves(root)
        fval = value(box)
        for i = 1:N
            bb = boxbounds(find_parent_with_splitdim(box, i), lower[i], upper[i])
            bb[2]-bb[1] < minwidth[i] && continue
            insert!(mes[i], width(box, i, x0, lower, upper), box=>extrapolate ? fval+qdelta(box, i) : fval)
        end
    end
    return mes
end

function trimschedule!(mes::Vector{<:MELink}, box::Box, splitdim, x0, lower, upper)
    for child in box.children
        fval = child.parent.fvalues[child.parent_cindex]
        for i = 1:ndims(box)
            trim!(mes[i], width(child, i, x0, lower, upper), child=>fval)
        end
    end
    return mes
end

function sweep!(root::Box, f, x0, splits, lower, upper; extrapolate::Bool = true, minwidth=zeros(eltype(x0), ndims(root)))
    mes = minimum_edges(root, x0, lower, upper, minwidth; extrapolate=extrapolate)
    sweep!(root, mes, f, x0, splits, lower, upper; minwidth=minwidth)
end
function sweep!(root::Box, mes::Vector{<:MELink}, f, x0, splits, lower, upper; minwidth=zeros(eltype(x0), ndims(root)))
    xtmp = similar(x0)
    flag = similar(x0, Bool)
    nsplits = similar(x0, Int)
    nleaves0 = count(x->true, leaves(root))
    nprocessed = 0
    visited = Set{typeof(root)}()
    dimorder = sortperm(length.(mes))  # process the dimensions with the shortest queues first
    for i in dimorder
        me = mes[i]
        while !isempty(me)
            item = popfirst!(me)
            box = item.l
            bb = boxbounds(find_parent_with_splitdim(box, i), lower[i], upper[i])
            bb[2]-bb[1] < minwidth[i] && continue
            position!(xtmp, flag, box)
            default_position!(xtmp, flag, x0)
            count_splits!(nsplits, box)
            if nsplits[i] > 0 && any(iszero, nsplits)
                # println("discarding ", box, " along dimension ", i)
                continue
            end
            nprocessed += 1
            empty!(visited)
            autosplit!(box, mes, f, x0, xtmp, i, splits, lower, upper, minwidth, visited)
        end
    end
    # println(nprocessed, " processed, starting with ", nleaves0, " leaves and ending with ", count(x->true, leaves(root)))
    root
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
    analyze!(root, f, x0, splits, lower, upper; rtol=rtol, atol=atol, fvalue=fvalue, maxevals=maxevals, kwargs...)
    return root, x0
end

"""
    root = analyze!(root, f, x0, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)

Further refinement of `root`. See [`analyze`](@ref) for details.
"""
function analyze!(root::Box, f::Function, x0, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500, print_interval=typemax(Int), kwargs...)
    box = minimum(root)
    boxval = value(box)
    lastval = typemax(boxval)
    tol_counter = 0
    lastprint = 0
    len = lenold = length(leaves(root))
    if print_interval < typemax(Int)
        println("Initial minimum ($len evaluations): ", minimum(root))
    end
    extrapolate = true
    while boxval > fvalue && tol_counter <= ndims(box) && len < maxevals
        lastval = boxval
        sweep!(root, f, x0, splits, lower, upper; extrapolate=extrapolate, kwargs...)
        extrapolate = !extrapolate
        box = minimum(root)
        boxval = value(box)
        len = length(leaves(root))
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
    end
    root
end

"""
    xmin, fmin = minimize(f, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)

Return the position `xmin` and value `fmin` of the minimum of `f` over the specified domain.
See [`analyze`](@ref) for information about the input arguments.
"""
function minimize(f, splits, lower, upper; rtol=1e-3, atol=0.0, fvalue=-Inf, maxevals=2500)
    root, x0 = analyze(f, splits, lower, upper; rtol=rtol, atol=atol, fvalue=fvalue, maxevals=maxevals)
    box = minimum(root)
    return position(box, x0), value(box)
end
