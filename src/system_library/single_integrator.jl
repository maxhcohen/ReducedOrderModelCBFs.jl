struct SingleIntegrator <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

SingleIntegrator() = SingleIntegrator(2, 2, x -> zeros(2), x -> diagm(ones(2)))
SingleIntegrator(N) = SingleIntegrator(N, N, x -> N == 1 ? 0.0 : zeros(N), x -> N == 1 ? 1.0 : diagm(ones(N)))
SingleIntegrator(::Type{ControlAffineSystem}) = SingleIntegrator()