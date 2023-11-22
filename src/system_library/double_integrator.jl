struct DoubleIntegratorControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end

struct DoubleIntegratorRobotic <: RoboticSystem
    n::Int
    m::Int
    M::Function
    H::Function
    B::Function
    C::Function
end

DoubleIntegrator(::Type{ControlAffineSystem}, N=2) = DoubleIntegratorControlAffine(N)
DoubleIntegrator(::Type{RoboticSystem}, N=2) = DoubleIntegratorRobotic(N)

function DoubleIntegratorControlAffine(N::Int)
    n = 2*N
    m = N
    f(x) = N == 1 ? [x[2], 0.0] : vcat(x[N+1:2*N], zeros(N))
    g(x) = N == 1 ? [0.0, 1.0] : vcat(zeros(N,N), diagm(ones(N)))

    return DoubleIntegratorControlAffine(n, m, f, g)
end

function DoubleIntegratorRobotic(N::Int)
    n = N
    m = N
    M(q) = n == 1 ? 1.0 : diagm(ones(N))
    C(q,v) = n == 1 ? 0.0 : zeros(N, N)
    G(q) = n == 1 ? 0.0 : zeros(N)
    H(q,v) = n == 1 ? 0.0 :  zeros(N)
    B(q) = diagm(ones(N))
    

    return DoubleIntegratorRobotic(n, m, M, H, B, C)
end

