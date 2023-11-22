struct CircularObstacle
    x
    r
end

function ControlBarrierFunction(Σ::SingleIntegrator, O::CircularObstacle, α::Function)
    h(x) = norm(x - O.x) - O.r^2

    return ControlBarrierFunction(h, α)
end
ControlBarrierFunction(Σ::SingleIntegrator, O::CircularObstacle) = ControlBarrierFunction(Σ, O, r -> r)

function ControlBarrierFunction(Σ::Unicycle, O::CircularObstacle, δ, α::Function)
    d(x) = norm(x[1:2] - O.x)
    θ(x) = atan(O.x[2] - x[2]. O.x[1] - x[1])
    h(x) = d(x) - O.r^2 - δ*cos(x[3] - θ(x))

    return ControlBarrierFunction(h, α)
end
ControlBarrierFunction(Σ::Unicycle, O::CircularObstacle) = ControlBarrierFunction(Σ, O, r -> r)