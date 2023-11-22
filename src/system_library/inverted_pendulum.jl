struct InvertedPendulumControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct InvertedPendulumRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    H::Function
    B::Function
end

function InvertedPendulum(
    ::Type{ControlAffineSystem};
    mass=1.0, 
    length=1.0, 
    damping=0.0,
)
    return InvertedPendulumControlAffine(mass=mass, length=length, damping=damping)
end

function InvertedPendulum(
    ::Type{RoboticSystem}; 
    mass=1.0, 
    length=1.0, 
    damping=0.0,
    ) 
    return InvertedPendulumRobotic(mass=mass, length=length, damping=damping)
end

function InvertedPendulumRobotic(; mass=1.0, length=1.0, damping=0.0)
    n = 1
    m = 1
    grav = 9.8
    M(q) = mass*length^2
    H(q, v) = mass*grav*length*sin(q) + damping*v
    B(q) = 1.0

    return InvertedPendulumRobotic(n, m, M, H, B)
end

function InvertedPendulumControlAffine(; mass=1.0, length=1.0, damping=0.0)
    Σ = InvertedPendulumRobotic(mass=mass, length=length, damping=damping)
    n, m, f, g = to_control_affine(Σ)

    return InvertedPendulumControlAffine(n, m, f, g)
end

