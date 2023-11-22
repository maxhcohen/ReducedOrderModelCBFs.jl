struct ControlBarrierFunction
    h::Function
    α::Function
end

(cbf::ControlBarrierFunction)(x) = cbf.h(x)

# Standard CBF constructorss
ControlBarrierFunction(h::Function) = ControlBarrierFunction(h, r -> r)
ControlBarrierFunction(h::Function, α::Float64) = ControlBarrierFunction(h, r -> α*r)

# CBF for robotic systems
function ControlBarrierFunction(Σ::RoboticSystem, h0::Function, k0::Function, α::Function; μ=1.0)
    # CBF for RoboticSystem
    h(q, v) = h0(q) - 0.5/μ * (v - k0(q))'*Σ.M(q)*(v - k0(q))

    # Add method for control affine system
    h(x) = Σ.n == 1 ? h(x[1], x[2]) : h(x[1:Σ.n], x[Σ.n+1:end])

    return ControlBarrierFunction(h, α)
end
