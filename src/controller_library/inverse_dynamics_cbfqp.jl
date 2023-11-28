"""
    InverseDynamicsCBFQP <: FeedbackController

WIP. Computes safe accelerations for a double integrator, which are translated into torques via inverse dynamics.
"""
struct InverseDynamicsCBFQP <: FeedbackController
    get_input::Function
end

(k::InverseDynamicsCBFQP)(q,q̇) = k.get_input(q, q̇)

"""
    InverseDynamicsCBFQP(Σ::RoboticSystem, cbf::ControlBarrierFunction, k0::SmoothSafetyFilter; kwargs)

Construct an inverse dynamics CBF-QP based controller.

# Arguments
- `Σ::RoboticSystem`: system under consideration
- `cbf::ControlBarrierFunction`: CBF for a reduced-order model
- `k0::SmoothSafetyFilter`: smooth safety filter for reduced-order model

# Keyword Arguments
- `μ=1.0`: parameter used to define CBF for full-order dynamics via backstepping
- `Kp=1.0`: proportional gain for tracking reduced-order model
- `Kd=1.0`: derivative gain for tracking reduced-order model
- `relaxed=false`: do we want to relax CBF condition on double integrator?
"""
function InverseDynamicsCBFQP(
    Σ::RoboticSystem, 
    cbf0::ControlBarrierFunction, 
    k0::SmoothSafetyFilter; 
    μ=1.0,
    Kp=1.0,
    Kd=1.0,
    relaxed=false,
)
    # Build control affine double integrator of same dimension as robotic system
    Σ1 = DoubleIntegrator(ControlAffineSystem, Σ.n)

    # Get robot dynamics
    M = Σ.M
    H = Σ.H
    B = Σ.B

    # Make CBF for double integrator
    h0 = cbf0.h
    e(x) = Σ.n == 1 ? x[2] - k0(x[1]) : x[Σ.n+1:end] - k0(x[1:Σ.n])
    V(x) = Σ.n == 1 ? e(x)'*M(x[1])*e(x) : e(x)'*M(x[1:Σ.n])*e(x)
    h(x) = Σ.n == 1 ? h0(x[1]) - (0.5/μ)*V(x) : h0(x[1:Σ.n]) - (0.5/μ)*V(x)
    Lfh, Lgh = get_lie_derivatives(Σ1, h)
    α = cbf0.α

    # Differentiate reduced-order controller
    dk0(q) = Σ.n == 1 ? ForwardDiff.derivative(k0, q) : ForwardDiff.jacobian(k0, q)
    q̈d(q,q̇) = dk0(q)*q̇

    # Now make input function
    function get_input(q, q̇)
        # Full-order state
        x = vcat(q, q̇)

        # Instantiate QP
        model = Model(OSQP.Optimizer)
        set_silent(model)

        # Make decision variables - acceleration q̈ and control input u
        q̈ = Σ.n == 1 ? @variable(model, q̈) : @variable(model, q̈[1:Σ.n])
        u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
        δ = @variable(model, δ)

        # Add CBF constraint for double integrator
        if relaxed
            @constraint(model, δ ≥ 0.0)
            @constraint(model, cbf_con, Lfh(x) + Lgh(x)*q̈ ≥ -α(h(x)) - δ)
        else
            @constraint(model, cbf_con, Lfh(x) + Lgh(x)*q̈ ≥ -α(h(x)))
        end

        # Add dynamics constraint to link q̈ and u
        @constraint(model, dyn_con, M(q)*q̈ + H(q,q̇) .== B(q)*u)

        # Add objective function - minimize control input
        @objective(model, Min, 0.5*u'u + 0.5*(q̈ + Kp*(q̇-k0(q)) + Kd*q̈d(q,q̇))'*(q̈ + Kp*(q̇-k0(q)) + Kd*q̈d(q,q̇)))

        # Solve
        optimize!(model)

        # Return control input
        return Σ.m == 1 ? value(u) : value.(u)
    end

    return InverseDynamicsCBFQP(get_input)
end