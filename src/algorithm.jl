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
# Returns `box, used_quasinewton`.
function autosplit!(box::Box{T}, mes::Vector{<:MELink}, f::CountedFunction, x0, xtmp, splitdim, xsplitdefaults, lower, upper, minwidth, visited::Set) where T
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
        return box, false
    end
    bb = boxbounds(p)
    bb[2]-bb[1] >= minwidth[splitdim] || return (box, false)
    N = ndims(box)
    nqn = ((N+1)*(N+2))÷2  # number of points needed for quasi-Newton approach
    if f.evals > 3*nqn     # make sure there is some excess
        Q, xbase, c = build_quadratic_model(box, x0)
        if Q.nzrows[] == size(Q.coefs, 1)
            g, B = solve(Q)
            success = quasinewton!(box, mes, B, g, c, f, x0, splitdim, lower, upper)
            return box, success
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
            return box, false
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
        (!success || nbr ∈ visited) && return (box, false) # don't get into a cycle
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
    return box, false
end

## Use regression to compute the best-fit quadratic model (with a dense Hessian)

# This implements the logic of "transposed" Gaussian elimination
# when adding a new point. We do it in transposed form because incoming
# points tend to build a lower-triangular matrix: the tree structure
# adds dimensions one at a time, so Δx will be all-zeros for trailing dimensions
# (once permuted into the order specified by dimpiv).
function Base.insert!(Q::QmIGE{T}, Δx, Δf, splitdim::Integer) where T
    myapprox(x, y, rtol) = isequal(x, y) | (abs(x-y) < rtol*(abs(x) + abs(y)))
    rtol = sqrt(eps(T))
    coefs, rhs, rowtmp = Q.coefs, Q.rhs, Q.rowtmp
    dimpiv, rowzero, ndims_old = Q.dimpiv, Q.rowzero, Q.ndims[]
    ndims = setrow!(rowtmp, dimpiv, Δx, ndims_old, splitdim)
    if ndims > ndims_old
        # Append to coefs. We always insert at the end of the
        # block corresponding to splitdim, so that we're performing
        # elimination on incoming points. That gives us a chance to
        # test for degeneracy before inserting them.
        i = ndims + (ndims*(ndims+1))÷2 # #gcoefs + #Bcoefs so far
        while i > 0 && rowtmp[i] == 0
            i -= 1
        end
        i == 0 && return Q
        for j = 1:i
            coefs[i,j] = rowtmp[j]
        end
        rowzero[i] = false
        rhs[i] = Δf
        Q.ndims[] = ndims
        Q.nzrows[] += 1
        return Q
    end
    # We've seen all the nonzero dimensions in Δx previously
    # Use elimination to determine whether it provides novel information
    lastnz = ndims
    while lastnz > 1 && Δx[dimpiv[lastnz]] == 0
        lastnz -= 1
    end
    i = lastnz + (lastnz*(lastnz+1))÷2
    @inbounds while i > 0
        if rowtmp[i] == 0
            i -= 1
            continue
        end
        rowzero[i] && break
        if abs(rowtmp[i]) > abs(coefs[i,i])
            # Swap (for numeric stability)
            for j = 1:i
                rowtmp[j], coefs[i,j] = coefs[i,j], rowtmp[j]
            end
            Δf, rhs[i] = rhs[i], Δf
        end
        c = rowtmp[i]/coefs[i,i]
        for j = 1:i-1
            sub = c * coefs[i,j]
            # rowtmp[j] = rowtmp[j] ≈ sub ? zero(sub) : rowtmp[j] - sub
            rowtmp[j] = myapprox(rowtmp[j], sub, rtol) ? zero(sub) : rowtmp[j] - sub
        end
        rowtmp[i] = 0  # in case of roundoff error
        Δf -= c*rhs[i]
        i -= 1
    end
    if i == 0
        return Q
    end
    for j = 1:i
        coefs[i,j] = rowtmp[j]
    end
    rowzero[i] = false
    rhs[i] = Δf
    Q.nzrows[] += 1
    return Q
end

function setrow!(rowtmp, dimpiv, Δx, ndims, splitdim)
    fill!(rowtmp, 0)
    k = 0
    have_splitdim = false
    for j = 1:ndims
        d = dimpiv[j]
        have_splitdim |= splitdim == d
        # The g coefs are linear in x, the B coefs quadratic
        # The implied storage order here matches `solve` below
        rowtmp[k+=1] = Δxd = Δx[d]  # g coef
        for i = 1:j-1
            rowtmp[k+=1] = Δxd * Δx[dimpiv[i]] # B[d,dimpiv[i]] coefficient
        end
        rowtmp[k+=1] = (Δxd * Δxd)/2 # B[d,d] coefficient
    end
    if !have_splitdim
        # We've not seen this dimension previously, so make it the
        # next one in dimpiv.
        ndims += 1
        dimpiv[ndims] = splitdim
        rowtmp[k+=1] = Δxd = Δx[splitdim]
        for i = 1:ndims-1
            rowtmp[k+=1] = Δxd * Δx[dimpiv[i]]
        end
        rowtmp[k+=1] = (Δxd * Δxd)/2
    end
    return ndims
end

function solve(Q::QmIGE{T,N}) where {T,N}
    gB = LowerTriangular(Q.coefs) \ Q.rhs
    g, B = Vector{T}(uninitialized, N), Matrix{T}(uninitialized, N, N)
    k = 0
    dimpiv = Q.dimpiv
    for i = 1:N
        di = dimpiv[i]
        g[di] = gB[k+=1]
        for j = 1:i
            dj = dimpiv[j]
            B[di,dj] = B[dj,di] = gB[k+=1]
        end
    end
    return g, B
end

function build_quadratic_model(box::Box{T,N}, x0) where {T,N}
    Q = QmIGE{T,N}()
    c = value(box)
    xbase = position(box, x0)
    if !isleaf(box)
        for i = 1:3
            descend!(Q, box.children[i], x0, xbase, c)
        end
    end
    while Q.nzrows[] < length(Q.rhs) && !isroot(box)
        cindex = box.parent_cindex
        box = box.parent
        j = 1
        if j == cindex j+=1 end
        descend!(Q, box.children[j], x0, xbase, c)
        j += 1
        if j == cindex j+=1 end
        descend!(Q, box.children[j], x0, xbase, c)
    end
    return Q, xbase, c
end

function descend!(Q, box, x0, xbase, c, skip=false)
    Q.nzrows[] == length(Q.rhs) && return Q
    Δx = position(box, x0)
    thisx = isleaf(box) ? zero(eltype(Δx)) : Δx[box.splitdim]
    if !skip
        for i = 1:length(Δx)
            Δx[i] -= xbase[i]
        end
        insert!(Q, Δx, value(box)-c, box.parent.splitdim)
    end
    if !isleaf(box)
        iskip = findfirst(equalto(thisx), box.xvalues)
        for i = 1:3
            descend!(Q, box.children[i], x0, xbase, c, i==iskip)
        end
    end
    return Q
end

function quasinewton!(box::Box{T}, mes, B, g, c, f, x0, splitdim, lower, upper, itermax = 20) where T
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
    # Do a backtracking linesearch
    fbox = value(box)
    while iter < itermax
        iter += 1
        xtarget = x + α*Δx
        # Sadly, the following function evaluations get discarded. But recording them
        # while preserving the box structure would take ndims(box)-1 additional
        # evaluations, which is not worth it unless xtarget itself is an improvement.
        if f(xtarget) > fbox
            α /= 2
            continue
        end
        leaf = find_leaf_at(root, xtarget)

        # # If leaf or one of its ancestors has been targeted before from an "external" box,
        # # terminate. The only allowed re-targetings are from inside the narrowest box yet
        # # targeted. This prevents running many line searches that all point to the same minimum.
        # if leaf != box
        #     dims_targeted = qtargeted(leaf, x, lower, upper)
        #     if all(dims_targeted)
        #         println("targeted ", xtarget, " with ", α)
        #         iter == 1 && record_targeted!(mes, leaf, splitdim)
        #         return false
        #     end
        # end

        # For consistency with the tree structure, we can change only one coordinate of
        # xleaf per split. Cycle through all coordinates, picking the one at each stage
        # that yields the smallest value as predicted by the quadratic model among the
        # not-yet-selected coordinates.
        dims_targeted = falses(ndims(leaf))
        q(x) = (x'*B*x)/2 + g'*x + c
        for j = 1:ndims(leaf)
            leaf.qtargeted = true
            xleaf = position(leaf, x0)
            xtest = copy(xleaf)
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
            if a == b || b == c
                # The box might be so narrow that there are not enough distinct floating-point
                # numbers that lie between the bounds.
                dims_targeted[imin] = true
                continue
            end
            split!(leaf, f, xleaf, imin, MVector3{T}(a, b, c), bb..., xleaf[imin], value(leaf))
            childindex = findfirst(equalto(xt), leaf.xvalues)
            leaf = leaf.children[childindex]
            dims_targeted[imin] = true
        end
        return true
    end
    return false
end

# Might want to add leaf to mes?
record_targeted!(mes, leaf, splitdim) = nothing

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
        fval = value(box)
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
    used_quasinewton = false
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
            _, qn = autosplit!(box, mes, f, x0, xtmp, i, splits, lower, upper, minwidth, visited)
            used_quasinewton |= qn
        end
    end
    # println(nprocessed, " processed, starting with ", nleaves0, " leaves and ending with ", count(x->true, leaves(root)))
    root, used_quasinewton
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
    fc = CountedFunction(f)
    box = minimum(root)
    boxval = value(box)
    lastval = typemax(boxval)
    tol_counter = 0
    lastprint = 0
    baseline_evals = len = lenold = length(leaves(root))
    if print_interval < typemax(Int)
        println("Initial minimum ($len evaluations): ", minimum(root))
    end
    extrapolate = true
    # The quasi-Newton step can reduce the function value so significantly, insist on
    # using it at least once.
    used_quasinewton = false
    while boxval > fvalue && (tol_counter <= ndims(box) || !used_quasinewton) && len < maxevals
        lastval = boxval
        _, qn = sweep!(root, fc, x0, splits, lower, upper; extrapolate=extrapolate, kwargs...)
        used_quasinewton |= qn
        extrapolate = !extrapolate
        box = minimum(root)
        boxval = value(box)
        len = baseline_evals + fc.evals
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
