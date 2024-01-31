struct PlanarQuadrotorControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct PlanarQuadrotorRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    C::Function
    G::Function
    H::Function
    B::Function
end

PlanarQuadrotor(::Type{ControlAffineSystem}) = PlanarQuadrotorControlAffine()
PlanarQuadrotor(::Type{RoboticSystem}) = PlanarQuadrotorRobotic()

function PlanarQuadrotorRobotic()
    # Quadrotor parameters
    g = 9.81
    m = 1.0
    J = 0.25

    # Inertia matrix
    M(q) = diagm([m, m, J])

    # Coriolis matrix
    C(q, q̇) = zeros(3,3)

    # Gravity vector
    G(q) = [0.0, m*g, 0.0]

    # Dynamics bias
    H(q, q̇) = C(q,q̇)*q̇ + G(q)

    # Actuation matrix
    B(q) = [sin(q[3]) 0.0; cos(q[3]) 0.0; 0.0 -1.0]


    return PlanarQuadrotorRobotic(3, 2, M, C, G, H, B)
end

function PlanarQuadrotorControlAffine()
    Σ = PlanarQuadrotorRobotic()
    n, m, f, g = to_control_affine(Σ)

    return PlanarQuadrotorControlAffine(n, m, f, g)
end

# Add extra stuff to stop from penetrating ground
function simulate(Σ::PlanarQuadrotorRobotic, q0, v0, T)
    x0 = ArrayPartition(q0, v0)
    n = Σ.n
    function odefun(dx, x, p, t)
        q = n == 1 ? x.x[1][1] : x.x[1]
        v = n == 1 ? x.x[2][1] : x.x[2]
        dx.x[1] .= n == 1 ? [v] : v
        dx.x[2] .= n == 1 ? [dynamics(Σ, q, v)] : dynamics(Σ, q, v)
    end

    condition(u, t, integrator) = u.x[1][2]
    affect!(integrator) = integrator.u.x[1][2] = 0.0
    cb = ContinuousCallback(condition, affect!)

    return solve(ODEProblem(odefun, x0, (0,T)), callback=cb)
end