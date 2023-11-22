struct CartPoleControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct CartPoleRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    C::Function
    H::Function
    G::Function
    B::Function
end

function CartPole(::Type{RoboticSystem}; mc=1.0, mp=0.2, l=0.5)
    return CartPoleRobotic(mc=mc, mp=mp, l=l)
end

function CartPole(::Type{ControlAffineSystem}; mc=1.0, mp=0.2, l=0.5)
    return CartPoleControlAffine(mc=mc, mp=mp, l=l)
end

function CartPoleRobotic(; mc=1.0, mp=0.2, l=0.5)
    # Inertia matrix
    M11(q) = mc + mp
    M12(q) = mp*l*cos(q[2])
    M21(q) = M12(q)
    M22(q) = mp*l^2
    M(q) = [M11(q) M12(q); M21(q) M22(q)]

    # Coriolis matrix
    C(q, q̇) = [0.0 -mp*l*q̇[2]*sin(q[2]); 0.0 0.0]

    # Gravitational effects
    g = 9.8
    G(q) = [0.0, mp*g*l*sin(q[2])]

    # Dynamics bias
    H(q, q̇) = C(q, q̇)*q̇ + G(q)

    # Actuation matrix
    B(q) = [1.0, 0.0]

    return CartPoleRobotic(2, 1, M, C, H, G, B)
end

function CartPoleControlAffine(; mc=1.0, mp=0.2, l=0.5)
    Σ = CartPoleRobotic(mc=mc, mp=mp, l=l)
    n, m, f, g = to_control_affine(Σ)

    return CartPoleControlAffine(n, m, f, g)
end