struct TunableCBFQP <: QPController
    get_input::Function
end

(k::TunableCBFQP)(x) = k.get_input(x)

function TunableCBFQP(Σ::ControlAffineSystem, cbf::ControlBarrierFunction, αd::Float64, kd::Function; relax=false)
    # System dimensions
    m = Σ.m

    # Lie derivative functions
    h = cbf.h
    Lfh, Lgh = get_lie_derivatives(Σ, cbf)

    function get_input(x)
        # Make optimizer
        model = Model(OSQP.Optimizer)
        set_silent(model)

        # Make decision variable for control input
        u = m == 1 ? @variable(model, u) : @variable(model, u[1:m])

        # Make decision variable for class K function coefficient
        α = @variable(model, α)

        # Set cost function
        @objective(model, Min, 0.5*(u - kd(x))'*(u - kd(x)) + 0.5*(α - αd)^2)

        # Set CBF constraint
        @constraint(model, Lfh(x) + Lgh(x)*u ≥ -α*h(x))

        # Make sure coefficients are positive
        if !relax
            @constraint(model, α ≥ 0.0)
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

    return TunableCBFQP(get_input)
end

function TunableCBFQP(
    Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, αds::Vector{Float64}, kd::Function;
    relax=false
)
    # System dimensions
    m = Σ.m
    N = length(cbfs)

    # Lie derivative functions
    hs = [cbf.h for cbf in cbfs]
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

        # Make decision variable for class K function coefficient
        αs = @variable(model, αs[1:N])

        # Set cost function
        @objective(model, Min, 0.5*(u - kd(x))'*(u - kd(x)) + sum([0.5*(α - αd)^2 for (α,αd) in zip(αs,αds)]))

        # Set CBF constraints
        for (Lfh, Lgh, h, α) in zip(Lfhs, Lghs, hs, αs)
            @constraint(model, Lfh(x) + Lgh(x)*u ≥ -α*h(x))
        end

        # Make sure coefficients are positive
        if !relax
            for α in αs
                @constraint(model, α ≥ 0.0)
            end
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

    return TunableCBFQP(get_input)
end

function TunableCBFQP(Σr::RoboticSystem, cbfs::ControlBarrierFunction, αds::Float64, kd::Function; relax=false)
    # Convert robotic system to control affine
    n, m, f, g = to_control_affine(Σr)
    Σ = CustomControlAffineSystem(n, m, f, g)
    return TunableCBFQP(Σ, cbfs, αds, kd, relax=relax)
end

function TunableCBFQP(Σr::RoboticSystem, cbfs::Vector{ControlBarrierFunction}, αds::Vector{Float64}, kd::Function; relax=false)
    # Convert robotic system to control affine
    n, m, f, g = to_control_affine(Σr)
    Σ = CustomControlAffineSystem(n, m, f, g)
    return TunableCBFQP(Σ, cbfs, αds, kd, relax=relax)
end