"""
    ReluSafetyFilter <: SafetyFilter

CBF-based safety filter resulting from solving a CBF-QP with a single constraint.
"""
struct ReluSafetyFilter <: SafetyFilter
    get_input::Function
end
(k::ReluSafetyFilter)(x) = k.get_input(x)
(k::ReluSafetyFilter)(q, q̇) = k.get_input(q, q̇)

"""
    ReluSafetyFilter(Σ::ControlAffineSystem, h::Function, α::Function, kd::Function; ε=0.0)
    ReluSafetyFilter(Σ::ControlAffineSystem, h::Function, α::Function; ε=0.0)
    ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction, kd::Function; ε=0.0)
    ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction; ε=0.0)
    ReluSafetyFilter(Σr::RoboticSystem, hr::Function, α::Function, kd::Function; ε=0.0)

Construct a ReluSafetyFilter for a control affine or robotic system.
"""
function ReluSafetyFilter(Σ::ControlAffineSystem, h::Function, α::Function, kd::Function; ε=0.0)
    # Get control affine dynamics
    n = Σ.n
    f = Σ.f
    g = Σ.g

    # Compute Lie derivatives
    ∇h(x) = n == 1 ? ForwardDiff.derivative(h, x) : ForwardDiff.gradient(h, x) 
    Lfh(x) = ∇h(x)'f(x)
    Lgh(x) = ∇h(x)'g(x)
    a(x) = ε == 0.0 ? Lfh(x) + Lgh(x)*kd(x) + α(h(x)) : Lfh(x) + Lgh(x)*kd(x) + α(h(x)) - (1/ε)*norm(Lgh(x))^2
    b(x) = norm(Lgh(x))^2

    # Solution of CBF-QP
    k(x) = kd(x) + λRelu(a(x), b(x))*Lgh(x)'

    return ReluSafetyFilter(k)
end

function ReluSafetyFilter(Σ::ControlAffineSystem, h::Function, α::Function; ε=0.0)
    # If not passed a desired controller, set it equal to zero and then construct controller
    kd(x) = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    return ReluSafetyFilter(Σ, h, α, kd, ε=ε)
end

function ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction, kd::Function; ε=0.0)
    # Pull out CBF functions
    h = cbf.h
    α = cbf.α

    return ReluSafetyFilter(Σ, h, α, kd, ε=ε)
end

function ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction; ε=0.0)
    # Make zeroing desired controller
    kd(x) = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    return ReluSafetyFilter(Σ, cbf, kd, ε=ε)
end

function ReluSafetyFilter(Σr::RoboticSystem, hr::Function, α::Function, kdr::Function; ε=0.0)
    # Dimension of configuration
    N = configuration_dim(Σr)

    # Convert robotic system to control affine system
    Σ = CustomControlAffineSystem(Σr)

    # Make sure CLF is compatible with control affine representation
    h(x) = N == 1 ? hr(x[1], x[2]) : hr(x[1:N], x[N+1:end])

    # Make sure desired controller us compatible with control affine representation
    kd(x) = N == 1 ? kdr(x[1], x[2]) : kdr(x[1:N], x[N+1:end])

    # Now make CBF controller for control affine system
    k = ReluSafetyFilter(Σ, h, α, kd, ε=ε)

    # Add method to handle configuration and velocity input
    get_input(q, q̇) = k.get_input(vcat(q, q̇))

    # Return new CBF controller with modified input function

    return ReluSafetyFilter(get_input)
end

# Some helper functions
relu(x) = max(0,x)
λRelu(a, b) = b == 0.0 ? 0.0 : relu(-a/sqrt(b))/sqrt(b)