"""
    CLFMinNorm <: FeedbackController

Pointwise Min Norm CLF controller.
"""
struct CLFMinNorm <: FeedbackController
    get_input::Function
end
(k::CLFMinNorm)(x) = k.get_input(x)
(k::CLFMinNorm)(q, q̇) = k.get_input(q, q̇)

"""
    CLFMinNorm(Σ::ControlAffineSystem, V::Function, α::Function)
    CLFMinNorm(Σr::RoboticSystem, Vr::Function, αr::Function)

Construct a CLF controller based on the pointwise min norm formula.
"""
function CLFMinNorm(Σ::ControlAffineSystem, V::Function, α::Function)
    # Extract control affine dynamics
    f = Σ.f
    g = Σ.g

    # Compute Lie derivatives
    ∇V(x) = Σ.n == 1 ? ForwardDiff.derivative(V, x) : ForwardDiff.gradient(V, x)
    LfV(x) = ∇V(x)'f(x)
    LgV(x) = ∇V(x)'g(x)
    a(x) = LfV(x) + α(x)
    b(x) = norm(LgV(x))^2

    # Pointwise min-norm controller
    k(x) = -λRelu(-a(x), b(x))*LgV(x)'

    return CLFMinNorm(k)
end

function CLFMinNorm(Σr::RoboticSystem, Vr::Function, αr::Function)
    # Dimension of configuration
    N = configuration_dim(Σr)

    # Convert robotic system to control affine system
    Σ = CustomControlAffineSystem(Σr)

    # Make sure CLF is compatible with control affine representation
    V(x) = N == 1 ? Vr(x[1], x[2]) : Vr(x[1:N], x[N+1:end])
    α(x) = N == 1 ? αr(x[1], x[2]) : αr(x[1:N], x[N+1:end])

    # Now make CLF controller for control affine system
    k = CLFMinNorm(Σ, V, α)

    # Add method to handle configuration and velocity input
    get_input(q, q̇) = k.get_input(vcat(q, q̇))

    # Reutn new CLF controller with modified input function

    return CLFMinNorm(get_input)
end