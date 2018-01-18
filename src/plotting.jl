using Plots, Colors, PerceptualColourMaps, Compat
using QuadDIRECT
using QuadDIRECT: Box

get_finite(x, default) = isfinite(x) ? x : default
outerbounds(b1, b2) = (min(b1[1], b2[1]), max(b1[2], b2[2]))
outerbounds(b1, b2::Real) = (min(b1[1], b2), max(b1[2], b2))
outerbounds_finite(b1, b2) = (min(b1[1], get_finite(b2[1], b1[1])), max(b1[2], get_finite(b2[2], b1[2])))
widenbounds(b) = (Δx = 0.1*(b[2]-b[1]); return (b[1]-Δx, b[2]+Δx))

function plotbounds(root::Box{T,2}, x0, lower, upper) where T
    xbounds, ybounds = (Inf, -Inf), (Inf, -Inf)
    bb = Vector{Tuple{T,T}}(uninitialized, 2)
    x = [NaN, NaN]
    flag = [false, false]
    for box in QuadDIRECT.visit_leaves(root)
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

function plotboxes(root::Box{T,2}, x0, clim, cs::AbstractVector{<:Colorant}, lower, upper, xbounds::Tuple{Any,Any}, ybounds::Tuple{Any,Any}) where T
    plt = plot(clims=clim, xlims=xbounds, ylims=ybounds)
    for bx in QuadDIRECT.visit_leaves(root)
        plt = plotrect!(plt, bx, x0, clim, cs, lower, upper, xbounds, ybounds)
    end
    plt
end
function plotboxes(root::Box{T,2}, x0, lower, upper, xbounds::Tuple{Any,Any}, ybounds::Tuple{Any,Any}) where T
    clim = extrema(root)
    cs = cmap("L4")
    plotboxes(root, x0, clim, cs, lower, upper, xbounds, ybounds)
end
plotboxes(root::Box{T,2}, x0, lower, upper) where T = plotboxes(root, x0, lower, upper, plotbounds(root, x0, lower, upper)...)

function plotrect!(plt, box::Box{T,2}, x0, clim, cs, lower, upper, xbounds, ybounds) where T
    @assert(QuadDIRECT.isleaf(box))
    x = QuadDIRECT.position(box, x0)
    bb = [(lower[i], upper[i]) for i = 1:2]
    QuadDIRECT.boxbounds!(bb, box)
    bbx, bby = bb
    rectx = [bbx[1], bbx[2], bbx[2], bbx[1], bbx[1]]
    recty = [bby[1], bby[1], bby[2], bby[2], bby[1]]
    fval = box.parent.fvalues[box.parent_cindex]
    cnorm = (clamp(fval, clim...) - clim[1])/(clim[2] - clim[1])
    col = cs[round(Int, (length(cs)-1)*cnorm) + 1]
    plot!(plt, clamp.(rectx, xbounds...), clamp.(recty, ybounds...), fill=(0,0.5,col), linecolor=:black, linealpha=0)
    plot!(plt, rectx, recty, linecolor=:black)
    scatter!(plt, [x[1]], [x[2]], markershape=:circle, markercolor=:black)
    plt
end
