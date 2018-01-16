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
