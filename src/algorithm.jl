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
        if xsplit[l] == xcur
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

function autosplit!(box::Box{T}, f, xtmp, splitdim, xsplitdefault, lower::Real, upper::Real) where T
    if !isleaf(box)
        # This box already got split along a different dimension
        box = find_smallest_child_leaf(box)
    end
    p = find_parent_with_splitdim(box, splitdim)
    isroot(p) && p.splitdim != splitdim && return split!(box, f, xtmp, splitdim, xsplitdefault, lower, upper)
    bb = boxbounds(p, lower, upper)
    xcur = p.parent.xvalues[p.parent_cindex]
    fcur = box.parent.fvalues[p.parent_cindex]
    # Do a quadratic fit
    xp, fp = p.parent.xvalues, p.parent.fvalues
    Δx = 2*max(xp[2]-xp[1], xp[3]-xp[2]) # in case the bounds are infinite, grow the gap
    xvert, fval, qcoef = qfit(xp[1]=>fp[1], xp[2]=>fp[2], xp[3]=>fp[3])
    if qcoef > 0 && xvert != xcur
        # xvert is a minimum
        if bb[1] <= xvert <= bb[2]
            # xvert is in the box
            a, b = lohi(xcur, xvert)
            imin, imax = biggest_interval(bb[1], a, b, bb[2])
            if isinf(imin)
                a, b, c = a-Δx, a, b
            elseif isinf(imax)
                a, b, c = a, b, b+Δx
            else
                a, b, c = lohi(a, b, (imin+imax)/2)
            end
            return split!(box, f, xtmp, splitdim, MVector3{T}(a, b, c), bb..., xcur, fcur)
        end
        #xvert is not in the box. TODO? should we use xvert somehow? (split a different box?)
    end
    # Trisect
    l, r = bb
    if isinf(l)
        l = r - 3*Δx  # gap is usually 1/3 of the box size, so undo this at the level of "box size"
    elseif isinf(r)
        r = l + 3*Δx
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
    for box in visit_leaves(root)
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
    for (i, me) in enumerate(mes)
        for item in me
            box = me.l
            position!(xtmp, flag, box)
            default_position!(xtmp, flag, x0)
            autosplit!(me.l, f, xtmp, i, splits[i], lower[i], upper[i])
        end
    end
    root
end
