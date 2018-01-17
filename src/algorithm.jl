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
    x0 = Tx[x[2] for x in xsplits]
    _init(f, x0, xsplits, lower, upper)
end

@noinline function _init(f, x0, xsplits, lower, upper)
    T = promote_type(eltype(x0), eltype(lower), eltype(upper))
    root = box = Box{T}()
    n = length(x0)
    xtmp = copy(x0)
    for i = 1:n
        xtmp = ipcopy!(xtmp, x0)
        box = split!(box, f, xtmp, i, xsplits[i], lower, upper)
        x0 = replacecoordinate!(x0, i, box.parent.xvalues[box.parent_cindex])
    end
    box
end

function split!(box::Box{T}, f, xtmp, splitdim, xsplit, lower, upper) where T
    # Evaluate f along the splits, keeping track of the best
    fsplit = Vector{T}(uninitialized, 3)
    fmin, idxmin = convert(T, Inf), 0
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
    add_children!(box, splitdim, copy(xsplit), fsplit, lower, upper)
    return box.children[idxmin]
end
