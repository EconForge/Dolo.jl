function old(x::Int)
    print("Old behavior")
end

function new(x::Int)
    print("New behavior")
end

@deprecate old(x::Int) new(x::Int) false

old(2)