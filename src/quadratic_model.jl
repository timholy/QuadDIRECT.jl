## Use regression to compute the best-fit quadratic model (with a dense Hessian)

const scaledthresh = 100 # FIXME hardcoded

"""
    g, B = solve(Q::QmIGE)

Solve a quadratic model `Q` for the gradient `g` and Hessian `B`. See
[`build_quadratic_model`](@ref) and [`complete_quadratic_model!`](@ref) to build `Q`.
"""
function solve(Q::QmIGE{T,N}) where {T,N}
    # Importantly, this defines the storage order in Q
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

# Uses available points to set up the regression problem
"""
    Q = build_quadratic_model(box, x0, scale, thresh)

Start building a quadratic model of the function, using points already stored
in the tree that contains `box`. `scale` (see [`boxscale`](@ref)) is used in part to
ensure numeric stability, and will prevent inclusion of points that are too distant along
any coordinate axis.

`Q` is not guaranteed to be complete; the tree may not contain enough points, or enough
"independent" points, to determine all the parameters of `Q`. See also [`complete_quadratic_model!`](@ref).
"""
function build_quadratic_model(box::Box{T,N}, x0, scale, thresh) where {T,N}
    Q = QmIGE{T,N}()
    c = value(box)
    xbase = position(box, x0)
    xtmp = similar(xbase)
    flag = similar(xtmp, Bool)
    nanfill = fill(T(NaN), N)
    bbs = boxbounds(box, nanfill, nanfill)
    succeeded = true
    @assert(isleaf(box))
    if !isleaf(box)
        for i = 1:3
            succeeded &= descend!(Q, c, box.children[i], x0, xbase, xtmp, flag, bbs, scale, thresh)
        end
    end
    # To increase the numeric stability of the result, we do more boxes than
    # nominally needed (see the `swap` portion of `insert!`)
    need_extra = true
    while succeeded && (Q.nzrows[] < length(Q.rhs) || need_extra) && !isroot(box)
        if Q.nzrows[] == length(Q.rhs)
            need_extra = minimum(abs.(diag(Q.coefs))) < thresh  # because of scaling this doesn't need units
        end
        cindex = box.parent_cindex
        box = box.parent
        succeeded &= all(isfinite, box.fvalues) && !box.qnconverged
        succeeded || break
        j = 1
        if j == cindex j+=1 end
        succeeded &= descend!(Q, c, box.children[j], x0, xbase, xtmp, flag, bbs, scale, thresh)
        succeeded || break
        j += 1
        if j == cindex j+=1 end
        succeeded &= descend!(Q, c, box.children[j], x0, xbase, xtmp, flag, bbs, scale, thresh)
    end
    return Q, xbase, c, succeeded
end

# Evaluates the function at new points until the quadratic model is fully specified
"""
    complete_quadratic_model!(Q::QmIGE, c, box::Box, f::Function, x0, xbase, scale, splits, lower::AbstractVector, upper::AbstractVector, thresh)

Split boxes in the tree containing `box` until all the parameters of the quadratic model can
be determined. It is allowed to fail if one or more boxes are too small to be split further,
so it is advisable to check that all diagonal coefficients in `Q` exceed `thresh` before
proceeding.

See also [`solve`](@ref).
"""
function complete_quadratic_model!(Q::QmIGE, c, box::Box{T,N}, f::Function, x0, xbase, scale, splits, lower::AbstractVector, upper::AbstractVector, thresh) where {T,N}
    @assert(all(x->x>0, Q.dimpiv))
    xtmp  = similar(x0)
    Δx    = similar(x0)
    Δxtmp = similar(x0)
    flag  = similar(x0, Bool)
    dimpiv = Q.dimpiv
    rowoffset = 0
    # Sanity check
    i = 0
    for idim = 1:N
        for irow = 0:idim
            if isassigned(Q.rowbox, i+=1)
                @assert(isleaf(Q.rowbox[i]))
                position!(xtmp, flag, Q.rowbox[i], x0)
                for j = idim+1:N
                    sd = Q.dimpiv[j]
                    @assert(xtmp[sd]==xbase[sd])
                end
            end
        end
    end
    # Process block-by-block, each block corresponding to a "new" splitting dimension
    for idim = 1:N
        all_are_set = true
        for j = 1:idim+1
            all_are_set &= abs(Q.coefs[j+rowoffset,j+rowoffset]) >= thresh
        end
        if all_are_set
            rowoffset += idim + 1
            continue
        end
        rowb = row = rowoffset + idim + 1  # within this dimension's block, go backwards
        irow = idim
        # @show Q.rowbox[rowoffset+1:row]
        # If all boxes that differed along this dimension were out-of-scale, we
        # need to split along that dimension.
        if !isassigned(Q.rowbox, rowb)
            splitdim = Q.dimpiv[idim]
            position!(xtmp, flag, box, x0)
            xcur, s = xbase[splitdim], scale[splitdim]
            p = find_parent_with_splitdim(box, splitdim)
            bb = boxbounds(p)
            @assert xtmp[splitdim] == xcur
            l = max(bb[1], xcur - s/2)
            u = min(bb[2], xcur + s/2)
            if l >= xcur
                l, xcur = xcur, (xcur + u)/2
            elseif xcur >= u
                u, xcur = xcur, (l + xcur)/2
            end
            split!(box, f, xtmp, splitdim, MVector3{T}(l, xcur, u), lower[splitdim], upper[splitdim], xcur, value(box))
            xcur = xbase[splitdim]
            j = 3
            j1 = 1;    if box.xvalues[j1] == xcur j = j1; j1 += 1; end
            j2 = j1+1; if box.xvalues[j2] == xcur j = j2; j2 += 1; end
            nbrbox = box.fvalues[j1] < box.fvalues[j2] ? box.children[j1] : box.children[j2]
            Q.rowbox[rowb] = nbrbox
            box = box.children[j]
        end
        while irow >= 0
            # @show row
            # if row == rowb
            #     rng = rowoffset+1:rowoffset+idim+1
            #     display(Q.coefs[rng,rng])
            #     @show thresh
            # end
            Qd = abs(Q.coefs[row,row])
            if Qd < thresh
                # @show rowb irow
                rbox = Q.rowbox[rowb]
                @assert isleaf(rbox)
                position!(xtmp, flag, rbox, x0)
                Δx .= xtmp .- xbase
                if irow > 0
                    # When irow > 0, splitting along the corresponding dimension will fill
                    # in this missing information
                    splitdim = dimpiv[irow]
                else
                    # root = get_root(box)
                    # splitprint_colored(root, box)
                    # println()
                    # splitprint_colored(root, rbox)
                    # println()
                    # @show Δx dimpiv
                    # irow == 0 corresponds to the coefficient of the gradient
                    if Δx[dimpiv[idim]] == 0
                        # There's no displacement in the gradient coordinate, so obviously
                        # we have to split it
                        splitdim = dimpiv[idim]
                    else
                        # Here it's not so obvious which coordinate we need to split.
                        # To avoid the risk of failure, we check the result of the gaussian
                        # elimination and choose our splitting dimension as the one that
                        # has largest (remaining) diagonal value.
                        Δxtmp[:] .= Δx ./ scale
                        maxdiag, splitdim = zero(T), 0
                        rng = rowoffset+1:rowoffset+idim+1
                        Q.rowzero[row] = true
                        for ii = 1:idim
                            # The precise magnitude of the step doesn't matter for these "probe" trials
                            sd = dimpiv[ii]
                            xold = Δxtmp[sd]
                            Δxtmp[sd] = Δxtmp[sd] == 0 ? 1 : 2*Δxtmp[sd] # ensure a non-zero entry no matter what
                            # @show Δxtmp
                            setrow!(Q.rowtmp, dimpiv, Δxtmp, N, sd)
                            Δxtmp[sd] = xold
                            # @show ii sd
                            rtmp, _ = incremental_gaussian_elimination!(Q.rowtmp, Q.coefs, Q.rowzero, rng; canswap=false, debug=false)
                            thisdiag = abs(Q.rowtmp[row])
                            if thisdiag > maxdiag
                                maxdiag = thisdiag
                                splitdim = sd
                            end
                        end
                        if splitdim == 0
                            rng = rowoffset+1:rowoffset+idim+1
                            display(Q.coefs[rng,rng])
                            error("no direction works")
                        end
                        # @show maxdiag splitdim
                    end
                end
                xcur = xtmp[splitdim]
                p = find_parent_with_splitdim(rbox, splitdim)
                # Check that the box is big enough to split along this dimension
                if !isroot(p)
                    bb = boxbounds(p)
                    if !(bb[2] - bb[1] > epswidth(bb))
                        return Q  # fail (FIXME?)
                    end
                end
                insrows, rboxn = split_insert!(Q, rbox, p, row:rowoffset+idim+1, f, xbase, c, xtmp, splitdim, scale, splits, lower, upper)
                if Qd == 0 && (maximum(insrows) != row || Q.coefs[row,row] == 0)
                    # Debugging
                    println("Something went wrong")
                    @show numevals(f) insrows row irow rowb rowoffset idim Q.dimpiv
                    @show rbox.parent.splitdim box.parent.splitdim
                    @show splitdim boxbounds(rbox, lower, upper)
                    rng = rowoffset+1:rowoffset+idim+1
                    display(Q.coefs[rng,rng])
                    @show xcur rbox.xvalues (rbox.xvalues .- xcur)./scale[splitdim]
                end
                Qd == 0 && @assert(maximum(insrows) == row)
                if row == rowb
                    Q.rowbox[row] = rboxn  # ensure it's a leaf
                end
            end
            if isassigned(Q.rowbox, row) && isleaf(Q.rowbox[row])
                rowb = row
            end
            row -=1; irow -= 1
        end
        rowoffset += idim + 1
    end
    return Q
end

function split_insert!(Q, rbox::Box{T}, p, rng, f, xbase, c, xtmp, splitdim, scale, splits, lower, upper) where T
    debug = false
    debug && @show rng
    row = first(rng)
    xcur = xtmp[splitdim]
    if isroot(p)
        xsplit = MVector3{T}(splits[splitdim])
        @assert(xsplit[1] < xsplit[2] < xsplit[3])
        lwr, upr = lower[splitdim], upper[splitdim]
    else
        bb = boxbounds(p)
        lwr, upr = bb[1], bb[2]
        xcur, scalesd = xtmp[splitdim], scale[splitdim]
        lo, hi = max(lwr, 5*lwr/6 + xcur/6, xcur - scalesd), min(upr, xcur/6 + 5*upr/6, xcur + scalesd)
        if lo + sqrt(eps(T))*scalesd >= xcur
            lo, xcur = xcur, (xcur+hi)/2
        elseif xcur + sqrt(eps(T))*scalesd >= hi
            xcur, hi = (xcur+lo)/2, xcur
        end
        xsplit = MVector3{T}(lo, xcur, hi)
        xcur = xtmp[splitdim]
        @assert(xsplit[1] < xsplit[2] < xsplit[3])
    end
    j1 = 1
    if xsplit[j1] == xcur j1 += 1 end
    j2 = j1+1
    if xsplit[j2] == xcur j2 += 1 end
    l, h = xsplit[j1], xsplit[j2]
    Qd = abs(Q.coefs[row, row])
    if Qd != 0
        # This is a case where the existing entry is below thresh. Before proceeding, check
        # whether the new one would be an improvement---rbox might be too small along
        # splitdim to raise the value.
        xtmp[splitdim] = l
        setrow!(Q.rowtmp, Q.dimpiv, (xtmp-xbase)./scale, ndims(rbox), splitdim)
        incremental_gaussian_elimination!(Q.rowtmp, Q.coefs, Q.rowzero, rng; canswap=false, debug=false)
        debug && @show abs(Q.rowtmp[row])
        if abs(Q.rowtmp[row]) < Qd
            xtmp[splitdim] = h
            setrow!(Q.rowtmp, Q.dimpiv, (xtmp-xbase)./scale, ndims(rbox), splitdim)
            incremental_gaussian_elimination!(Q.rowtmp, Q.coefs, Q.rowzero, rng; canswap=false, debug=false)
            debug && @show abs(Q.rowtmp[row])
            if abs(Q.rowtmp[row]) < Qd
                # Neither was an improvement, so keep what we've got
                debug && println("skipping")
                return (row,-1), rbox
            end
        end
    end
    # Split rbox and insert into Q
    val = value(rbox)
    fsplit = MVector3(val,val,val)
    xtmp[splitdim] = l
    fsplit[j1] = f(xtmp)
    xtmp[splitdim] = h
    fsplit[j2] = f(xtmp)
    add_children!(rbox, splitdim, xsplit, fsplit, lwr, upr)
    if l != xbase[splitdim]
        xtmp[splitdim] = l
        debug && @show (xtmp.-xbase)./scale
        insrow1 = insert!(Q, (xtmp.-xbase)./scale, rbox.fvalues[j1]-c, splitdim, rbox.children[j1]; canswap=false, debug=debug)
    else
        debug && println("skipping l")
        insrow1 = 0
    end
    if h != xbase[splitdim]
        xtmp[splitdim] = h
        debug && @show (xtmp.-xbase)./scale
        insrow2 = insert!(Q, (xtmp.-xbase)./scale, rbox.fvalues[j2]-c, splitdim, rbox.children[j2]; canswap=false, debug=debug)
    else
        debug && println("skipping h")
        insrow2 = 0
    end
    # Return the leaf corresponding to the original
    xtmp[splitdim] = xcur
    j = 1; if j == j1 j += 1 end; if j == j2 j += 1 end
    return ((insrow1, insrow2), rbox.children[j])
end

function descend!(Q::QmIGE, c, box, x0, xbase, Δx, flag, bbs, scale, thresh, skip::Bool=false)
    box.qnconverged && return false  # don't build quadratic models that include points tagged as converged minima
    # Is this box so far away that we should prune this branch?
    boxbounds!(bbs, flag, box)
    for i = 1:ndims(box)
        xb, bb = xbase[i], bbs[i]
        min(bb[1] - xb, xb - bb[2]) > scaledthresh*scale[i] && return true
    end
    position!(Δx, flag, box, x0)
    thisx = isleaf(box) ? zero(eltype(Δx)) : Δx[box.splitdim]
    if !skip
        # Add this point to the model
        isgood = true
        leaf = find_leaf_at(box, Δx)
        for i = 1:length(Δx)
            Δx[i] = (Δx[i] - xbase[i])/scale[i]
            # To prevent distant points from numerically overwhelming the local structure of the box,
            # put in a cutoff. The concern is that large values here will dominate
            # eigenvalues that come from local data.
            isgood &= abs(Δx[i]) < scaledthresh
        end
        if isgood
            insert!(Q, Δx, value(box)-c, box.parent.splitdim, leaf)
        else
            # we have to save space for this dimension, since they are added
            # in splitdim order
            sd = box.parent.splitdim
            if sd ∉ Q.dimpiv
                Q.ndims[] += 1
                Q.dimpiv[Q.ndims[]] = sd
            end
        end
    end
    # Recursively add all children
    if !isleaf(box)
        iskip = findfirst(equalto(thisx), box.xvalues)
        for i = 1:3
            descend!(Q, c, box.children[i], x0, xbase, Δx, flag, bbs, scale, thresh, i==iskip) || return false
        end
    end
    return true
end

function Base.insert!(Q::QmIGE, Δx, Δf, splitdim::Integer, box::Box; canswap::Bool=true, debug::Bool=false)
    coefs, rhs, rowtmp = Q.coefs, Q.rhs, Q.rowtmp
    dimpiv, rowzero, ndims_old = Q.dimpiv, Q.rowzero, Q.ndims[]
    ndims = setrow!(rowtmp, dimpiv, Δx, ndims_old, splitdim)
    if ndims > ndims_old
        # Append to coefs. We always insert at the end of the
        # block corresponding to splitdim, so that we're performing
        # elimination on incoming points. That gives us a chance to
        # test for degeneracy before inserting them.
        i = ndims + (ndims*(ndims+1))÷2 # #gcoefs + #Bcoefs so far
        i == 0 && return i
        for j = 1:i
            coefs[i,j] = rowtmp[j]
        end
        rowzero[i] = false
        Q.rowbox[i] = box
        rhs[i] = Δf
        Q.ndims[] = ndims
        Q.nzrows[] += 1
        return i
    end
    ndims == 0 && return 0
    # We've seen all the nonzero dimensions in Δx previously
    # Use elimination to determine whether it provides novel information
    blockstart = ndims + ((ndims-1)*ndims)÷2    # first row for this splitting dimension
    blockend = blockstart + ndims               # last row for this splitting dimension
    finalrow, newrhs, box = incremental_gaussian_elimination!(rowtmp, coefs, rowzero, blockstart:blockend, Δf, rhs, box, Q.rowbox; canswap=canswap, debug=debug)
    if finalrow > 0
        rowzero[finalrow] = false
        for j = 1:finalrow
            coefs[finalrow,j] = rowtmp[j]
        end
        store!(rhs, newrhs, finalrow)
        store!(Q.rowbox, box, finalrow)
        Q.nzrows[] += 1
    end
    return finalrow
end

# This implements the logic of "flipped" Gaussian elimination (the zeros block on the
# upper triangle) when adding a new point. We do it in flipped form because incoming
# points tend to build a lower-triangular matrix: the tree structure
# adds dimensions one at a time, so Δx will be all-zeros for trailing dimensions
# (once permuted into the order specified by dimpiv).
function incremental_gaussian_elimination!(newrow::AbstractVector{T}, coefs::AbstractMatrix, rowzero::AbstractVector{Bool}, rng::AbstractUnitRange, newrhs=nothing, rhs=nothing, meta=nothing, metastore=nothing; canswap::Bool=true, debug::Bool=false) where T
    # @inline myapprox(x, y, rtol) = isequal(x, y) | (abs(x-y) < rtol*(abs(x) + abs(y)))
    @inline myapprox(x, y, rtol) = (x == y) | (abs(x-y) < rtol*(abs(x) + abs(y)))
    rtol = 1000*eps(T)

    debug && @show rng
    debug && println("input: ", newrow[rng])
    i = last(rng)
    @assert(first(rng) > 0)
    @inbounds while i >= first(rng)
        if newrow[i] == 0
            i -= 1
            continue
        end
        rowzero[i] && break
        if canswap && abs(newrow[i]) > abs(coefs[i,i])
            # Swap (for numeric stability)
            debug && println("swap at ", i, "($(newrow[i]) vs $(coefs[i,i])")
            for j = 1:i
                newrow[j], coefs[i,j] = coefs[i,j], newrow[j]
            end
            newrhs = swap!(rhs, newrhs, i)
            meta   = swap!(metastore, meta, i)
        end
        c = newrow[i]/coefs[i,i]
        @simd for j = 1:i-1
            sub = c * coefs[i,j] # coefsp[j,i]
            newrow[j] = ifelse(myapprox(newrow[j], sub, rtol), zero(sub), newrow[j] - sub)
        end
        newrow[i] = 0  # rather than subtracting, just set it (to avoid roundoff error)
        debug && @show i newrow[rng]
        newrhs = subtract(rhs, c, i, newrhs)
        i -= 1
    end
    if i < first(rng)
        i = 0
    end
    return i, newrhs, meta
end

function swap!(store, x, i)
    x, store[i] = store[i], x
    return x
end
swap!(::Any, ::Void, i) = nothing

store!(store, x, i) = store[i] = x
store!(::Any, ::Void, i) = nothing

subtract(rhs, c, i, newrhs) = newrhs - c*rhs[i]
subtract(rhs, c, i, ::Void) = nothing


function setrow!(rowtmp, dimpiv, Δx, ndims, splitdim)
    fill!(rowtmp, 0)
    k = 0
    have_splitdim = false
    for j = 1:ndims
        d = dimpiv[j]
        have_splitdim |= splitdim == d
        # The g coefs are linear in x, the B coefs quadratic
        # The implied storage order here matches `solve` below
        Δxd = Δx[d]
        if Δxd == 0
            k += j+1
            continue
        end
        rowtmp[k+=1] = Δxd  # g coef
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
    # Back up to last nonzero dimension
    while ndims > 0 && Δx[dimpiv[ndims]] == 0
        ndims -= 1
    end
    return ndims
end
