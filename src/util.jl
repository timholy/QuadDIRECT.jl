"""
   xvert, fvert, qcoef = qfit(xm=>fm, x0=>f0, xp=>fp)

Given three points `xm < x0 < xp ` and three corresponding
values `fm`, `f0`, and `fp`, fit a quadratic. Returns the position `xvert` of the vertex,
the quadratic's value `fvert` at `xvert`, and the coefficient `qcoef` of the quadratic term.
`xvert` is a minimum if `qcoef > 0`.

Note that if the three points lie on a line, `qcoef == 0` and both `xvert` and `fvert` will
be infinite.
"""
function qfit(xfm, xf0, xfp)
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    @assert(xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp))
    cm = fm/((xm-x0)*(xm-xp))  # coefficients of Lagrange polynomial
    c0 = f0/((x0-xm)*(x0-xp))
    cp = fp/((xp-xm)*(xp-x0))
    qcoef = cm+c0+cp
    qvalue(x) = cm*(x-x0)*(x-xp) + c0*(x-xm)*(x-xp) + cp*(x-xm)*(x-x0)
    if fm == f0 == fp
        return x0, f0, zero(qcoef)
    end
    xvert = (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*qcoef)
    return xvert, qvalue(xvert), qcoef
end

## Minimum Edge List utilities
Base.empty!(mel::MELink) = (mel.next = mel; return mel)

function Base.insert!(mel::MELink, w, lf::Pair)
    l, f = lf
    prev, next = mel, mel.next
    while prev != next && w > next.w
        if f < next.f
            if next != next.next
                # Drop the next item from the list
                next = next.next
                prev.next = next
            else
                # Drop the last item from the list
                prev.next = prev
                next = prev
            end
        else
            prev = next
            next = next.next
        end
    end
    if prev == next
        # we're at the end of the list
        prev.next = typeof(mel)(l, w, f)
    else
        if f < next.f
            prev.next = typeof(mel)(l, w, f, next)
        end
    end
    return mel
end

Base.start(mel::MELink) = mel
Base.done(mel::MELink, state::MELink) = state == state.next
Base.next(mel::MELink, state::MELink) = (state.next, state.next)

function Base.show(io::IO, mel::MELink)
    print(io, "List(")
    next = mel.next
    while mel != next
        print(io, '(', next.w, ", ", next.l, "=>", next.f, "), ")
        mel = next
        next = next.next
    end
    print(io, ')')
end
