struct Unicycle <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

function Unicycle()
    f(x) = zeros(3)
    g1(x) = [cos(x[3]), sin(x[3]), 0.0]
    g2(x) = [0.0, 0.0, 1.0]
    g(x) = hcat(g1(x), g2(x))

    return Unicycle(3, 2, f, g)
end

Unicycle(::Type{ControlAffineSystem}) = Unicycle()