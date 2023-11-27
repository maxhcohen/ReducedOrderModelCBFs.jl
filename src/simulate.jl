"""
    simulate(Σ::ControlAffineSystem, x0, T)
    simulate(Σ::ControlAffineSystem, k, x0, T)
    simulate(Σ::ControlAffineSystem, k::FeedbackController, x0, T)
    simulate(Σ::RoboticSystem, q0, v0, T)
    simulate(Σ::RoboticSystem, k, q0, v0, T)
    simulate(Σ::RoboticSystem, k::FeedbackController, q0, v0, T)

Simulate a system from an initial condition under a specified controller for a given time.
"""
function simulate(Σ::ControlAffineSystem, x0, T)
    function odefun(dx, x, p, t)
        dx[1:Σ.n] .= dynamics(Σ, x)
    end

    return solve(ODEProblem(odefun, x0, (0,T)))
end

function simulate(Σ::ControlAffineSystem, k, x0, T)
    function odefun(dx, x, p, t)
        dx[1:Σ.n] .= dynamics(Σ, x, k(x,t))
    end

    return solve(ODEProblem(odefun, x0, (0,T)))
end

function simulate(Σ::ControlAffineSystem, k::FeedbackController, x0, T)
    return simulate(Σ, (x, t) -> k(x), x0, T)
end

function simulate(Σ::RoboticSystem, q0, v0, T)
    x0 = ArrayPartition(q0, v0)
    n = Σ.n
    function odefun(dx, x, p, t)
        q = n == 1 ? x.x[1][1] : x.x[1]
        v = n == 1 ? x.x[2][1] : x.x[2]
        dx.x[1] .= n == 1 ? [v] : v
        dx.x[2] .= n == 1 ? [dynamics(Σ, q, v)] : dynamics(Σ, q, v)
    end

    return solve(ODEProblem(odefun, x0, (0,T)))
end

function simulate(Σ::RoboticSystem, k, q0, v0, T)
    x0 = ArrayPartition(q0, v0)
    n = Σ.n
    function odefun(dx, x, p, t)
        q = n == 1 ? x.x[1][1] : x.x[1]
        v = n == 1 ? x.x[2][1] : x.x[2]
        dx.x[1] .= n == 1 ? [v] : v
        dx.x[2] .= n == 1 ? [dynamics(Σ, q, v, k(q, v, t))] : dynamics(Σ, q, v, k(q, v, t))
    end

    return solve(ODEProblem(odefun, x0, (0,T)))
end

function simulate(Σ::RoboticSystem, k::FeedbackController, q0, v0, T)
    return simulate(Σ, (q,v,t) -> k(q, v), q0, v0, T)
end