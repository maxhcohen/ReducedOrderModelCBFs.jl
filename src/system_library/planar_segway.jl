struct PlanarSegwayControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct PlanarSegwayRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    H::Function
    B::Function
end

PlanarSegway(::Type{ControlAffineSystem}) = PlanarSegwayControlAffine()
PlanarSegway(::Type{RoboticSystem}) = PlanarSegwayRobotic()

function PlanarSegwayRobotic()
    # Segway parameters
    g = 9.81
    R = 0.195
    M = 2*2.485
    Jc = 2*0.0559
    L = 0.169
    m = 44.798
    Jg = 3.836
    m0 = 52.710
    J0 = 5.108
    Km = 2*1.262
    bt = 2*1.225

    # Dynamics
    D(q) = [m0 m*L*cos(q[2]); m*L*cos(q[2]) J0]
    H(q, q̇) = [-m*L*sin(q[2])*q̇[2] + bt*(q̇[1] - R*q̇[2])/R, -m*g*L*sin(q[2]) - bt*(q̇[1] - R*q̇[2])]
    B(q) = [Km/R, -Km]

    return PlanarSegwayRobotic(2, 1, D, H, B)
end

function PlanarSegwayControlAffine()
    Σ = PlanarSegwayRobotic()
    n, m, f, g = to_control_affine(Σ)

    return PlanarSegwayControlAffine(n, m, f, g)
end