"""
    mel = MELink{Tx,Ty}(dummylabel)

Create a linked-list storing labels and x/y pairs that represent the "minimum edge" of a
set of values, where the stored x and y values are monotonic in both x and y. Insertion
of a new x/y pair expunges any entries with smaller x but larger y. Add new elements with
`insert!(mel, x, label=>y)`.

The first item in the list is a dummy item with `x=0`, so that the head of the list is constant.
"""
mutable struct MELink{Tl,Tw,Tf}
    l::Tl   # label
    w::Tw   # x-coordinate (box width)
    f::Tf   # y-coordinate (function value)
    next::MELink{Tl,Tw,Tf}

    function MELink{Tl,Tw,Tf}(l, w, f) where {Tw, Tf, Tl}
        mel = new{Tl,Tw,Tf}(l, w, f)
        mel.next = mel
        return mel
    end
    function MELink{Tl,Tw,Tf}(l, w, f, next) where {Tw, Tf, Tl}
        mel = new{Tl,Tw,Tf}(l, w, f, next)
        return mel
    end
end

MELink{Tw,Tf}(dummylabel::Tl) where {Tl,Tw,Tf} = MELink{Tl,Tw,Tf}(dummylabel, 0, -Inf)

Base.iteratorsize(::Type{<:MELink}) = Base.SizeUnknown()
