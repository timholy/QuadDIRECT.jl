using PyPlot, Colors, PerceptualColourMaps
using QuadDIRECT
using QuadDIRECT: Box

get_finite(x, default) = isfinite(x) ? x : default
outerbounds(b1, b2) = (min(b1[1], b2[1]), max(b1[2], b2[2]))
outerbounds(b1, b2::Real) = (min(b1[1], b2), max(b1[2], b2))
outerbounds_finite(b1, b2) = (min(b1[1], get_finite(b2[1], b1[1])), max(b1[2], get_finite(b2[2], b1[2])))
widenbounds(b) = (Δx = 0.1*(b[2]-b[1]); return (b[1]-Δx, b[2]+Δx))

function plotbounds(root::Box{T,2}, x0, lower, upper) where T
    xbounds, ybounds = (Inf, -Inf), (Inf, -Inf)
    bb = Vector{Tuple{T,T}}(undef, 2)
    x = [NaN, NaN]
    flag = [false, false]
    for box in leaves(root)
        for i = 1:2
            bb[i] = (lower[i], upper[i])
        end
        QuadDIRECT.boxbounds!(bb, box)
        xbounds = outerbounds_finite(xbounds, bb[1])
        ybounds = outerbounds_finite(ybounds, bb[2])
        x = QuadDIRECT.position!(x, flag, box)
        QuadDIRECT.default_position!(x, flag, x0)
        xbounds = outerbounds(xbounds, x[1])
        ybounds = outerbounds(ybounds, x[2])
    end
    widenbounds(xbounds), widenbounds(ybounds)
end

function plotboxes(root::Box{T,2}, x0, lower, upper, xbounds::Tuple{Any,Any}, ybounds::Tuple{Any,Any}; clim=extrema(root), cs::AbstractVector{<:Colorant}=cmap("RAINBOW3")) where T
    clf()
    for bx in leaves(root)
        plotrect(bx, x0, clim, cs, lower, upper, xbounds, ybounds)
    end
    xlim(xbounds)
    ylim(ybounds)
    ax = gca()
    # norm = PyPlot.matplotlib[:colors][:Normalize](vmin=clim[1], vmax=clim[2])
    # cb = PyPlot.matplotlib[:colorbar][:ColorbarBase](ax, cmap=ColorMap(cs), norm=norm)
    ax
end
plotboxes(root::Box{T,2}, x0, lower, upper; clim=extrema(root), cs::AbstractVector{<:Colorant}=cmap("RAINBOW3")) where T =
    plotboxes(root, x0, lower, upper, plotbounds(root, x0, lower, upper)...; clim=clim, cs=cs)

function plotrect(box::Box{T,2}, x0, clim, cs, lower, upper, xbounds, ybounds) where T
    @assert(QuadDIRECT.isleaf(box))
    x = QuadDIRECT.position(box, x0)
    bb = QuadDIRECT.boxbounds(box, lower, upper)
    bbx, bby = bb
    rectx = [bbx[1], bbx[2], bbx[2], bbx[1], bbx[1]]
    recty = [bby[1], bby[1], bby[2], bby[2], bby[1]]
    fval = box.parent.fvalues[box.parent_cindex]
    cnorm = (clamp(fval, clim...) - clim[1])/(clim[2] - clim[1])
    col = cs[round(Int, (length(cs)-1)*cnorm) + 1]
    hold(true)
    fill(clamp.(rectx, xbounds...), clamp.(recty, ybounds...), color=[red(col), green(col), blue(col)])
    plot(rectx, recty, "k")
    plot([x[1]], [x[2]], "k.")
    hold(false)
    nothing
end
