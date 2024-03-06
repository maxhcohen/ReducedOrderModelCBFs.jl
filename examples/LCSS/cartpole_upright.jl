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
θd = π
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

# Get safety filer
kCBF = ReluSafetyFilter(Σr, h, r -> α*r, (q,q̇) -> 0.0)

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
sol2 = simulate(Σr, q0, q̇0, T)

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
        ax_theme...,
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
        xmin=0,
        xmax=5,
        ymin=θd - θmax - 0.1,
        ymax=θd + θmax + 0.1,
        legend_style={draw="none", fill="none", legend_cell_align="left", font="\\footnotesize", at="{(0.0,0.75)}",anchor="north west"},
        height="1.7in",
        width="1.9in",
    },
    Plot({plt_theme..., color=colors[3]}, Coordinates(ts, sol.(ts, idxs=2))),
    LegendEntry("CBF"),
    Plot({plt_theme..., color=colors[4], dotted, opacity=0.6}, Coordinates(ts, sol2.(ts, idxs=2))),
    LegendEntry("No CBF"),
    Plot({plt_theme..., dashed}, Coordinates([0, T], [θd + θmax, θd + θmax])),
    Plot({plt_theme..., dashed}, Coordinates([0, T], [θd - θmax, θd - θmax])),
    [raw"\filldraw[gray, thick, opacity=0.4] (0,2.88) -- (10,2.88) -- (10,2) -- (0, 2) -- cycle;"],
    [raw"\filldraw[gray, thick, opacity=0.4] (0,3.4) -- (10,3.4) -- (10,4) -- (0, 4) -- cycle;"],
    TextNode(3.0, 3.35, raw"{\footnotesize$\theta_{\mathrm{d}} +\theta_{\max}$}"),
    TextNode(2.3, 2.95, raw"{\footnotesize$\theta_{\mathrm{d}} -\theta_{\max}$}"),
)
# pgfsave("cartpole_upright_theta_small.pdf", fig2)


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

# pgfsave("cartpole_down.pdf", fig3)

# Plot safe set?
contour_opts = @pgf {
    contour_lua={labels=false,levels=0.0},
    line_width="1.0pt",
}

x1 = θd-θmax-0.2:0.01:θd+θmax+0.2
x2 = -2.0:0.1:2.0
fig4 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$\theta$",
        ylabel=raw"$\dot{\theta}$",
        xmin=x1[1],
        xmax=x1[end],
        ymin=x2[1],
        ymax=x2[end],
        # major_tick_length="0.075cm",
        view="{0}{90}",
        "colormap/blackwhite",
        height="1.7in",
        width="1.9in",
    },
    Plot3(
        {
            contour_opts...,
        },
        Table(x1, x2, [h([0.0,x], [0.0,y]) for x in x1, y in x2]),
    ),   
    Plot({plt_theme..., dashed}, Coordinates([θd-θmax, θd-θmax], [-3, 3])),
    Plot({plt_theme..., dashed}, Coordinates([θd+θmax, θd+θmax], [-3, 3])),
    [raw"\filldraw[gray, thick, opacity=0.4] (2.88,-3) -- (2.88,3) -- (2,3) -- (2, -3) -- cycle;"],
    [raw"\filldraw[gray, thick, opacity=0.4] (3.4,-3) -- (3.4,3) -- (4,3) -- (4, -3) -- cycle;"],
    # [raw"\filldraw[gray, thick, opacity=0.4] (-5.5,-1) -- (3.5,-1) -- (3.5,-0.262) -- (-5.5, -0.262) -- cycle;"],
)
# pgfsave("theta_safe_set.pdf", fig4)