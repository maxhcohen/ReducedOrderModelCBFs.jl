struct DoublePendulumControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct DoublePendulumRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    H::Function
    B::Function
    C::Function
    G::Function
end

function DoublePendulum(::Type{RoboticSystem}; m1=2.0, m2=2.0, l1=1.0, l2=1.0)
    return DoublePendulumRobotic(m1=m1, m2=m2, l1=l1, l2=l2)
end

function DoublePendulum(::Type{ControlAffineSystem}; m1=2.0, m2=2.0, l1=1.0, l2=1.0)
    return DoublePendulumControlAffine(m1=m1, m2=m2, l1=l1, l2=l2)
end

function DoublePendulumRobotic(;m1=2.0, m2=2.0, l1=1.0, l2=1.0)
    # Mass matrix
    M11(q) = (m1 + m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(q[2])
    M12(q) = m2*l2^2 + m2*l1*l2*cos(q[2])
    M22(q) = m2*l2^2
    M(q) = [M11(q) M12(q); M12(q) M22(q)]

    # Coriolis matrix
    C11(q, q̇) = 0.0
    C12(q, q̇) = -m2*l1*l2*(2*q̇[1] + q̇[2])*sin(q[2])
    C21(q, q̇) = -0.5*m2*l1*l2*(2*q̇[1] + q̇[2])*sin(q[2])
    C22(q, q̇) = -0.5*m2*l1*l2*q̇[1]*sin(q[2])
    C(q, q̇) = [C11(q,q̇) C12(q, q̇); C21(q, q̇) C22(q, q̇)]

    # Gravitational affects
    g = 9.8
    G1(q) = (m1 + m2)*l1*sin(q[1]) + m2*l2*sin(q[1] + q[2])
    G2(q) = m2*l2*sin(q[1] + q[2])
    G(q) = g*[G1(q), G2(q)]

    # Dynamic bias term
    H(q,q̇) = C(q,q̇)*q̇ + G(q)

    # Actuation matrix
    B(q) = [1.0 0.0; 0.0 1.0]

    return DoublePendulumRobotic(2, 2, M, H, B, C, G)
end

function DoublePendulumControlAffine(; m1=2.0, m2=2.0, l1=1.0, l2=1.0)
    Σ = DoublePendulumRobotic(m1=m1, m2=m2, l1=l1, l2=l2)
    n, m, f, g = to_control_affine(Σ)

    return DoublePendulumControlAffine(n, m, f, g)
end