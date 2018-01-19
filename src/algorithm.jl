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
    for i = 1:n
        xtmp = ipcopy!(xtmp, xstar)
        box = split!(box, f, xtmp, i, xsplits[i], lower[i], upper[i])
        xstar = replacecoordinate!(xstar, i, box.parent.xvalues[box.parent_cindex])
    end
    box, xstar
end

function split!(box::Box{T}, f, xtmp, splitdim, xsplit, lower::Real, upper::Real, xcur=NaN, fcur=NaN) where T
    # Evaluate f along the splits, keeping track of the best
    fsplit = MVector3{T}(typemax(T), typemax(T), typemax(T))
    fmin, idxmin = typemax(T), 0
    for l = 1:3
        @assert(isfinite(xsplit[l]))
        if xsplit[l] == xcur && !isnan(fcur)
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
    return box.children[idxmin]
end

function autosplit!(box::Box{T}, f, x0, xtmp, splitdim, xsplitdefaults, lower, upper, oldbox=Set([box])) where T
    if !isleaf(box)
        # This box already got split along a different dimension
        box = find_smallest_child_leaf(box)
        xtmp = position(box, x0)
    end
    if splitdim == 0
        # If we entered this box from a neighbor, just choose an axis that has been least-split
        nsplits = count_splits(box)
        splitdim = indmin(nsplits)
        if nsplits[splitdim] > 0
            # If all dimensions have been split, prefer splitting unbounded dimensions
            bbs = boxbounds(box, lower, upper)
            if all(isfinite, bbs[splitdim])
                for (i, bb) in enumerate(bbs)
                    if !isfinite(bb[1]) || !isfinite(bb[2])
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
    isroot(p) && p.splitdim != splitdim && return split!(box, f, xtmp, splitdim, xsplitdefault, lwr, upr)
    bb = boxbounds(p)
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
            return add_children!(box, splitdim, MVector3{T}(xf1[1], xf2[1], xf3[1]),
                                 MVector3{T}(xf1[2], xf2[2], xf3[2]), lwr, upr)
        end
        # xvert is not in the box. Prepare to split the neighbor, but for this box
        # just bisect xcur and xnew
        if xvert < bb[1]
            xtmp, dir = replacecoordinate!(xtmp, splitdim, bb[1]), -1
        else
            xtmp, dir = replacecoordinate!(xtmp, splitdim, bb[2]), +1
        end
        nbr = find_leaf_at_edge(get_root(box), xtmp, splitdim, dir)
        xmid = (xcur+xnew)/2
        xtmp = replacecoordinate!(xtmp, splitdim, xmid)
        fmid = f(xtmp)
        xf1, xf2, xf3 = order_pairs(xcur=>fcur, xnew=>fnew, xmid=>fmid)
        if (length(oldbox) == 1 && first(oldbox) == box) || nbr ∈ oldbox
            add_children!(box, splitdim, MVector3{T}(xf1[1], xf2[1], xf3[1]),
                        MVector3{T}(xf1[2], xf2[2], xf3[2]), lwr, upr)
        end
        nbr ∈ oldbox && return box # don't get into a cycle
        return autosplit!(nbr, f, x0, position(nbr, x0), 0, xsplitdefaults, lower, upper, push!(oldbox, box))
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
    return split!(box, f, xtmp, splitdim, MVector3{T}(a, b, c), bb..., xcur, fcur)
end

# A dumb O(N) algorithm for building the minimum-edge structures
function minimum_edges(root::Box{T,N}, x0, lower, upper) where {T,N}
    mes = [MELink{T,T}(root) for i = 1:N]
    for box in leaves(root)
        fval = box.parent.fvalues[box.parent_cindex]
        for i = 1:N
            insert!(mes[i], width(box, i, x0, lower, upper), box=>fval)
        end
    end
    return mes
end

function sweep!(root::Box, f, x0, splits, lower, upper)
    mes = minimum_edges(root, x0, lower, upper)
    xtmp = similar(x0)
    flag = similar(x0, Bool)
    nsplits = similar(x0, Int)
    nleaves0 = count(x->true, leaves(root))
    nprocessed = 0
    for (i, me) in enumerate(mes)
        for item in me
            box = item.l
            position!(xtmp, flag, box)
            default_position!(xtmp, flag, x0)
            count_splits!(nsplits, box)
            if nsplits[i] > 0 && any(iszero, nsplits)
                println("discarding ", box, " along dimension ", i)
                continue
            end
            nprocessed += 1
            autosplit!(box, f, x0, xtmp, i, splits, lower, upper)
        end
    end
    println(nprocessed, " processed, starting with ", nleaves0, " leaves and ending with ", count(x->true, leaves(root)))
    root
end
