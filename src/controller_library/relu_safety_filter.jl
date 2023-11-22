struct ReluSafetyFilter <: SafetyFilter
    get_input::Function
end

(k::ReluSafetyFilter)(x) = k.get_input(x)

function ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction, kd::Function; ε=0.0)
    f = Σ.f
    g = Σ.g
    h = cbf.h
    α = cbf.α
    ∇h(x) = ForwardDiff.gradient(h, x)
    Lfh(x) = ∇h(x)'f(x)
    Lgh(x) = ∇h(x)'g(x)
    a(x) = ε == 0.0 ? Lfh(x) + Lgh(x)*kd(x) + α(h(x)) : Lfh(x) + Lgh(x)*kd(x) + α(h(x)) - (1/ε)*norm(Lgh(x))^2
    b(x) = norm(Lgh(x))^2
    k(x) = kd(x) + λRelu(a(x), b(x))*Lgh(x)'

    return ReluSafetyFilter(k)
end

function ReluSafetyFilter(Σ::ControlAffineSystem, cbf::ControlBarrierFunction; ε=0.0)
    kd(x) = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    return ReluSafetyFilter(Σ, cbf, kd, ε=ε)
end

# Some helper functions
relu(x) = max(0,x)
λRelu(a, b) = b == 0.0 ? 0.0 : relu(-a/sqrt(b))/sqrt(b)