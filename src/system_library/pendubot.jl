struct PendubotControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct PendubotRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    C::Function
    H::Function
    G::Function
    B::Function
end

function Pendubot(::Type{RoboticSystem})
    return PendubotRobotic()
end

function Pendubot(::Type{ControlAffineSystem})
    return PendubotControlAffine()
end

function PendubotRobotic()
    # Inertia matrix
    M11(q) = 2.0
    M12(q) = cos(q[1] - q[2])
    M21(q) = M12(q)
    M22(q) = 1.0
    M(q) = [M11(q) M12(q); M21(q) M22(q)]

    # Coriolis matrix
    C(q, q̇) = [0.0 sin(q[1] - q[2])*q̇[2]; -sin(q[1] - q[2])*q̇[1] 0.0]

    # Gravitational effects
    g = 9.8
    P(q) = 2*g*cos(q[1]) + g*cos(q[2])
    G(q) = ForwardDiff.gradient(P, q)

    # Dynamics bias
    H(q, q̇) = C(q, q̇)*q̇ + G(q)

    # Actuation matrix
    B(q) = [1.0, 0.0]

    return PendubotRobotic(2, 1, M, C, H, G, B)
end

function PendubotControlAffine()
    Σ = PendubotRobotic()
    n, m, f, g = to_control_affine(Σ)

    return PendubotControlAffine(n, m, f, g)
end