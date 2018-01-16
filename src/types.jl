mutable struct MELink{Tl,Tw,Tf}
    l::Tl
    w::Tw
    f::Tf
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
