using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σr = CartPole(RoboticSystem)

# System output 1 - cart position
y1(q) = q[1]
J1(q) = [1.0, 0.0]'

# System output 2 - pendulum angle
y2(q) = q[2]
J2(q) = [0.0, 1.0]'

# Wall location
pmax = 2.0

# Output constraint 1
h01(y) = pmax - y

# Desired pendulum angle and maximum deviation
θd = 0.0
θmax = π/12

# Output constraint
h02(y) = θmax^2 - (θd - y)^2

# Desired controller for full-order system: move to a desired location
xd = 3.0
Kp = 1.0
Kd = 2.0
kd(q, q̇) = -Kp*(q[1] - xd) - Kd*q̇[1] + 0.05*q̇[2]

# Desired controller for reduced-order system
kd0(y) = 0.0

# Reduced-order CBF parameters for first CBF
α = 5.0
∇h01(y) = ForwardDiff.derivative(h01, y)
Lfh01(y) = 0.0
Lgh01(y) = ∇h01(y)
a1(y) = Lfh01(y) + Lgh01(y)*kd0(y) + α*h01(y)
b1(y) = norm(Lgh01(y))^2

# Reduced-order CBF parameters for second CBF
α = 5.0
∇h02(y) = ForwardDiff.derivative(h02, y)
Lfh02(y) = 0.0
Lgh02(y) = ∇h02(y)
a2(y) = Lfh02(y) + Lgh02(y)*kd0(y) + α*h02(y)
b2(y) = norm(Lgh02(y))^2

# Smooth safety filters
σ = 0.05
k01(q) = λSoftplus(a1(y1(q)), b1(y1(q)), σ)*Lgh01(y1(q))'
k02(q) = λSoftplus(a2(y2(q)), b2(y2(q)), σ)*Lgh02(y2(q))'

# Now build CBF sfor full-order system
μ = 10.0
h1(q,q̇) = h01(y1(q)) - (0.5/μ)*norm(J1(q)*q̇ - k01(q))^2
h2(q,q̇) = h02(y2(q)) - (0.5/μ)*norm(J2(q)*q̇ - k02(q))^2
h1(x) = h1(x[1:2], x[3:4])
h2(x) = h2(x[1:2], x[3:4])
cbf1 = ControlBarrierFunction(h1)
cbf2 = ControlBarrierFunction(h2)

# Make tunable CBFQP 
kTCBF = TunableCBFQP(Σr, [cbf1, cbf2], [1.0, 1.0], x -> kd(x[1:2], x[3:4]))
k(q, q̇) = kTCBF([q;q̇])

# Initial conditions
p0 = -5.0
θ0 = 0.2
ṗ0 = 0.0
θ̇0 = 0.0
q0 = [p0, θ0]
q̇0 = [ṗ0, θ̇0]

# Sim parameters
T = 10.0
dt = 0.05
ts = 0.0:dt:T

# Simulate
sol = simulate(Σr, (q,q̇,t) -> k(q,q̇), q0, q̇0, T)

# Plot results
ax_theme = get_ieee_theme_medium()
plt_theme = get_plt_theme()
colors = get_colors()
fig1 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$x(t)$",
        ax_theme...,
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, sol.(ts, idxs=1))),
    Plot({plt_theme...,}, Coordinates([0,T], [pmax, pmax])),
)

fig2 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "thick"}, Coordinates([0, T], [θmax, θmax])),
    Plot({"smooth", "thick"}, Coordinates([0, T], -[θmax, θmax])),
)

fig3 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$u(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, k.(sol.(ts, idxs=1:2), sol.(ts, idxs=3:4)))),
)

# fig4 = @pgf Axis(
#     {
#         ax_theme...,
#         width="2.4in",
#         height="2.1in",
#         xlabel=raw"$x$",
#         ylabel=raw"$\theta$",
#         xmin=-5.5,
#     },
#     Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
#     Plot({plt_theme...,}, Coordinates([pmax, pmax], [-θmax, θmax])),
#     Plot({plt_theme...,}, Coordinates([-6, pmax], [θmax, θmax])),
#     Plot({plt_theme...,}, Coordinates([-6, pmax], -[θmax, θmax])),
# )

fig4 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$\theta$",
        xmin=-5.5,
        xmax=3.5,
        ymin=-1.0,
        ymax=1.0,
    },
    [raw"\filldraw[gray, thick, opacity=0.0, fill opacity=0.4] (2,1) -- (3.5,1) -- (3.5,-1) -- (2, -1) -- cycle;"],
    [raw"\filldraw[gray, thick, opacity=0.0, fill opacity=0.4] (-5.5,1) -- (2,1) -- (2,0.262) -- (-5.5, 0.262) -- cycle;"],
    [raw"\filldraw[gray, thick, opacity=0.0, fill opacity=0.4] (-5.5,-1) -- (2,-1) -- (2,-0.262) -- (-5.5, -0.262) -- cycle;"],
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    Plot({plt_theme...,}, Coordinates([pmax, pmax], [-θmax, θmax])),
    Plot({plt_theme...,}, Coordinates([-6, pmax], [θmax, θmax])),
    Plot({plt_theme...,}, Coordinates([-6, pmax], -[θmax, θmax])),
    TextNode(p0+0.5, θ0-0.1, raw"{$\mathbf{q}_0$}"),
    TextNode(-1, -0.5, raw"{$-\theta_{\max}$}"),
    TextNode(-1, 0.5, raw"{$\theta_{\max}$}"),
    TextNode(2.4, -0.5, raw"{$x_{\max}$}"),
)

# pgfsave("cartpole_combined.pdf", fig4)