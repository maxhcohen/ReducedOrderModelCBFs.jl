# Abstract system types
abstract type System end
abstract type ControlAffineSystem <: System end
abstract type TwoLayerControlAffineSystem <: ControlAffineSystem end
abstract type RoboticSystem <: ControlAffineSystem end

# Get state and control dimensions
state_dim(Σ::ControlAffineSystem) = Σ.n
state_dim(Σ::RoboticSystem) = 2*Σ.n
configuration_dim(Σ::RoboticSystem) = Σ.n
control_dim(Σ::System) = Σ.m
is_fully_actuated(Σ::RoboticSystem) = configuration_dim(Σ) == control_dim(Σ)

# Dynamics functions
dynamics(Σ::ControlAffineSystem, x) = Σ.f(x)
dynamics(Σ::ControlAffineSystem, x, u) = Σ.f(x) + Σ.g(x)*u
dynamics(Σ::RoboticSystem, q, v) = -Σ.M(q)\Σ.H(q,v)
dynamics(Σ::RoboticSystem, q, v, u) = Σ.M(q)\(Σ.B(q)*u - Σ.H(q,v)) 
inverse_dynamics(Σ::RoboticSystem, q, v, v̇) = Σ.M(q)*v̇ + Σ.H(q,v)

# Construct a control affine system from a robotic one
"""
    to_control_affine(Σ::RoboticSystem)

Convert a robotic system into control affine form.
"""
function to_control_affine(Σ::RoboticSystem)
    n = 2*Σ.n
    m = Σ.m
    function f(x)
        if Σ.n == 1
            q = x[1]
            v = x[2]
            return [v, -Σ.M(q)\Σ.H(q,v)]
        else
            q = x[1:Σ.n]
            v = x[Σ.n+1:end]
            return vcat(v, -Σ.M(q)\Σ.H(q,v))
        end
    end

    function g(x)
        if Σ.n == 1
            q = x[1]
            v = x[2]
            return [0.0, Σ.M(q)\Σ.B(q)]
        else
            q = x[1:Σ.n]
            v = x[Σ.n+1:end]
            if Σ.m == 1
                return vcat(zeros(Σ.n), Σ.M(q)\Σ.B(q))
            else
                return vcat(zeros(Σ.n, Σ.m), Σ.M(q)\Σ.B(q))
            end
        end
    end

    return n, m, f, g
end
