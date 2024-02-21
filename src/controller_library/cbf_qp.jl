struct CBFQP <: QPController
    get_input::Function
end

(k::CBFQP)(x) = k.get_input(x)

function CBFQP(Σ::ControlAffineSystem, cbf::ControlBarrierFunction, kd::Function)
    # System dimensions
    m = Σ.m

    # Lie derivative functions
    h = cbf.h
    α = cbf.α
    Lfh, Lgh = get_lie_derivatives(Σ, cbf)

    function get_input(x)
        # Make optimizer
        model = Model(OSQP.Optimizer)
        set_silent(model)

        # Make decision vriable
        u = m == 1 ? @variable(model, u) : @variable(model, u[1:m])

        # Set cost function
        @objective(model, Min, 0.5*(u - kd(x))'*(u - kd(x)))

        # Set CBF constraint
        @constraint(model, Lfh(x) + Lgh(x)*u ≥ -α(h(x)))

        # Solve QP
        optimize!(model)

        # Return input
        if m == 1
            return value(u)
        else
            return value.(u)
        end
    end

    return CBFQP(get_input)
end

function CBFQP(Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)
    # System dimensions
    m = Σ.m

    # Lie derivative functions
    hs = [cbf.h for cbf in cbfs]
    αs = [cbf.α for cbf in cbfs]
    Lfhs = []
    Lghs = []
    for cbf in cbfs
        Lfh, Lgh = get_lie_derivatives(Σ, cbf)
        push!(Lfhs, Lfh)
        push!(Lghs, Lgh)
    end

    function get_input(x)
        # Make optimizer
        model = Model(OSQP.Optimizer)
        set_silent(model)

        # Make decision vriable
        u = m == 1 ? @variable(model, u) : @variable(model, u[1:m])

        # Set cost function
        @objective(model, Min, 0.5*(u - kd(x))'*(u - kd(x)))

        # Set CBF constraints
        for (Lfh, Lgh, α, h) in zip(Lfhs, Lghs, αs, hs)
            @constraint(model, Lfh(x) + Lgh(x)*u ≥ -α(h(x)))
        end

        # Solve QP
        optimize!(model)

        # Return input
        if m == 1
            return value(u)
        else
            return value.(u)
        end
    end

    return CBFQP(get_input)
end

function CBFQP(Σ::ControlAffineSystem, cbf::ControlBarrierFunction)
    # System dimensions
    m = Σ.m

    kd(x) = m == 1 ? 0.0 : zeros(m)
    
    return CBFQP(Σ, cbf, kd)
end