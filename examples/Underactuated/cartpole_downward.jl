using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σr = CartPole(RoboticSystem)

# System output - pendulum angle
y(q) = q[2]
J(q) = [0.0, 1.0]'

# Desired pendulum angle and maximum deviation
θd = 0.0
θmax = π/12

# Output constraint
h0(y) = θmax^2 - (θd - y)^2

# Reduced-order model: 1D single integrator
Σ0 = CustomControlAffineSystem(1, 1, x -> 0.0, x -> 1.0)

# Desired controller for reduced-order system: keep angle at θd
K = 1.0
kd0(y) = -K*(y - θd)

# Reduced-order CBF parameters
α = 5.0
k0 = SmoothSafetyFilter(Σ0, h0, r -> α*r, kd0)

# Now build CBF for full-order system
μ = 10.0
h(q,q̇) = h0(y(q)) - (0.5/μ)*norm(J(q)*q̇ - k0(y(q)))^2

# Desired controller for full-order system: move to a desired location
xd = 3.0
Kp = 1.0
Kd = 2.0
Kθ = 10.0
kd(q, q̇) = -Kp*(q[1] - xd) - Kd*q̇[1] + 0.1*q̇[2] + Kθ*(q[2] - k0(y(q)))

# Get safety filer
kCBF = ReluSafetyFilter(Σr, h, r -> α*r, kd)

# Check out decoupling matrix
D = Σr.M
B = Σr.B
A(q) = J(q)*inv(D(q))*B(q)

# Initial conditions
p0 = -5.0
θ0 = θd - 0.1
ṗ0 = 0.0
θ̇0 = 0.0
q0 = [p0, θ0]
q̇0 = [ṗ0, θ̇0]

# Sim parameters
T = 10.0
dt = 0.05
ts = 0.0:dt:T

# Simulate
sol = simulate(Σr, kCBF, q0, q̇0, T)

# Plot results
ax_theme = get_ieee_theme_medium()
plt_theme = get_plt_theme()
colors = get_colors()
fig1 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$p(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=1))),
)

fig2 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "thick"}, Coordinates([0, T], [θd + θmax, θd + θmax])),
    Plot({"smooth", "thick"}, Coordinates([0, T], [θd - θmax, θd - θmax])),
)

# fig3 = @pgf Axis(
#     {
#         ax_theme...,
#         xlabel=raw"$x$",
#         ylabel=raw"$\theta$",
#         width="2.4in",
#         height="2.1in",
#         xmin=-5.5,
#     },
#     Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
#     Plot({plt_theme...}, Coordinates([pmax, pmax], [θd - θmax, θd + θmax])),
#     Plot({plt_theme...}, Coordinates([-6, pmax], [θd + θmax, θd + θmax])),
#     Plot({plt_theme...}, Coordinates([-6, pmax], [θd - θmax, θd - θmax])),
# )

fig3 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$\theta$",
        xmin=-5.5,
        xmax=3.5,
        ymin=-1.0,
        ymax=1.0,
    },
    [raw"\filldraw[gray, thick, opacity=0.4] (-5.5,1) -- (3.5,1) -- (3.5,0.262) -- (-5.5, 0.262) -- cycle;"],
    [raw"\filldraw[gray, thick, opacity=0.4] (-5.5,-1) -- (3.5,-1) -- (3.5,-0.262) -- (-5.5, -0.262) -- cycle;"],
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    # Plot({plt_theme..., dashed}, Coordinates([pmax, pmax], [-2, 2])),
    Plot({plt_theme...,}, Coordinates([-6, 4], [θd + θmax, θd + θmax])),
    Plot({plt_theme...,}, Coordinates([-6, 4], [θd - θmax, θd - θmax])),
    TextNode(p0+0.3, θ0+0.1, raw"{$\mathbf{q}_0$}"),
    TextNode(-1, -0.5, raw"{$-\theta_{\max}$}"),
    TextNode(-1, 0.5, raw"{$\theta_{\max}$}"),
    # TextNode(1.2, -0.7, raw"{$x_{\max}$}"),
)

pgfsave("cartpole_down.pdf", fig3)