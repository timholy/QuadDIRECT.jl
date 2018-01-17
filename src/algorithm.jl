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

function split!(box::Box{T}, f, xtmp, splitdim, xsplit, lower::Real, upper::Real) where T
    # Evaluate f along the splits, keeping track of the best
    fsplit = MVector3{T}(typemax(T), typemax(T), typemax(T))
    fmin, idxmin = typemax(T), 0
    for l = 1:3
        xtmp = replacecoordinate!(xtmp, splitdim, xsplit[l])
        ftmp = f(xtmp)
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
