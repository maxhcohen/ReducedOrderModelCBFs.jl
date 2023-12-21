using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX
using LaTeXStrings

# Load in system - cartpole
Σ = CartPole(RoboticSystem)

# Decompose dynamics to build CBF for underactuated system
D(q) = Σ.M(q)
D11(q) = D(q)[1,1]
D12(q) = D(q)[1,2]
D21(q) = D(q)[2,1]
D22(q) = D(q)[2,2]

H(q, q̇) = Σ.H(q,q̇)
H1(q, q̇) = Σ.H(q,q̇)[1]
H2(q, q̇) = Σ.H(q,q̇)[2]

B = Σ.B(zeros(2))
B1 = B[1]
B2 = B[2]

# Transformed system dynamics
D̄(q) = D12(q) - D11(q)*pinv(D21(q))*D22(q)
H̄(q, q̇) = H1(q, q̇) - D11(q)*pinv(D21(q))*H2(q,q̇)
B̄ = B1

# Now make CBF for system
θmax = π/6
h0(q) = θmax^2 - (q[2] - π)^2
∇h0(q) = ForwardDiff.gradient(h0, q)
a0(q) = h0(q)
b0(q) = norm(∇h0(q))^2
k0(q) = λSoftplus(a0(q), b0(q), 0.1)*∇h0(q)
h(q, q̇) = h0(q) - 0.5*(q̇[2] - k0(q)[2])'*D̄(q)'*D̄(q)*(q̇[2] - k0(q)[2])
h(x) = h(x[1:2], x[3:4])
∇h(x) = ForwardDiff.gradient(h, x)

# Now get control affine representation
Σ2 = CartPole(ControlAffineSystem)
f(x) = Σ2.f(x)
g(x) = Σ2.g(x)
Lfh(x) = ∇h(x)'f(x)
Lgh(x) = dot(∇h(x), g(x))

# Desired controller - track a angle reference
function θd(t)
    if t < 5
        return π - π/4
    else
        return π + π/4
    end
end
kd(x, t) = -50*(x[2] - θd(t)) - 30*x[4]

# Safety filter
a(x, t) = Lfh(x) + Lgh(x)*kd(x, t) + h(x)
b(x) = norm(Lgh(x))^2
kqp(x, t) = kd(x, t) + λRelu(a(x, t), b(x))*Lgh(x)'

# Initial conditions
p0 = 0.0
θ0 = π
ṗ0 = 0.0
θ̇0 = 0.0
q0 = [p0, θ0]
q̇0 = [ṗ0, θ̇0]
x0 = vcat(q0, q̇0)

# Simulation for sanity check
T = 10.0
ts = 0.0:0.01:T
sol = simulate(Σ2, kqp, x0, T)

# Plot states
ax_theme, plt_theme, colors = plot_defaults()
ax1 = @pgf Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        height="1.75in",
        width="1.75in",
        # xlabel=L"t",
        # ylabel=L"\theta(t)",
    }, 
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [π+θmax, π+θmax])),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [π-θmax, π-θmax])),
    Plot({plt_theme..., color=colors[2], opacity=0.5}, Coordinates([0, 5], [π-π/4, π-π/4])),
    Plot({plt_theme..., color=colors[2], opacity=0.5}, Coordinates([5, T], [π+π/4, π+π/4])),
)
# pgfsave("cartpole1.pdf", ax1)

# Plot CBF
ax2 = @pgf Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        width="1.75in",
        height="1.75in",
    }, 
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, h0.(sol.(ts, idxs=1:2)))),
    Plot({plt_theme..., "black", "dashed"}, Coordinates([0,T], [0, 0])),
)
# pgfsave("cartpole2.pdf", ax2)