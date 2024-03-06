using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σr = CartPole(RoboticSystem)

# System output - cart position
y(q) = q[1]
J(q) = [1.0, 0.0]'

# Wall location
pmax = 2.0

# Output constraint
h0(y) = pmax - y

# Desired pendulum angle and maximum deviation
θd = 0.0
θmax = π/12

# Desired controller for full-order system: move to a desired location
xd = 3.0
Kp = 1.0
Kd = 2.0
kd(q, q̇) = -Kp*(q[1] - xd) - Kd*q̇[1] + 0.05*q̇[2]

# Reduced-order model: 1D single integrator
Σ0 = CustomControlAffineSystem(1, 1, x -> 0.0, x -> 1.0)

# Desired controller for reduced-order system
kd0(y) = 0.0

# Smooth safety filter
α = 5.0
k0 = SmoothSafetyFilter(Σ0, h0, r -> α*r, kd0)

# Now build CBF for full-order system
μ = 10.0
h(q,q̇) = h0(y(q)) - (0.5/μ)*norm(J(q)*q̇ - k0(y(q)))^2

# Get safety filer
kCBF = ReluSafetyFilter(Σr, h, r -> α*r, kd)

# Check out decoupling matrix
D = Σr.M
B = Σr.B
A(q) = J(q)*inv(D(q))*B(q)

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
sol = simulate(Σr, kCBF, q0, q̇0, T)
sol2 = simulate(Σr, (q,q̇,t) -> kd(q,q̇), q0, q̇0, T)

# Plot results
# ax_theme = get_ieee_theme_medium()
ax_theme = get_ieee_theme_small()
plt_theme = get_plt_theme()
colors = get_colors()
fig1 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$t$",
        ylabel=raw"$x(t)$",
        xmin=0,
        xmax=T,
        ymax=4.0,
        legend_pos="south east",
        height="1.7in",
        width="1.9in",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, sol.(ts, idxs=1))),
    LegendEntry("CBF"),
    Plot({plt_theme..., color=colors[2], dotted, opacity=0.6,}, Coordinates(ts, sol2.(ts, idxs=1))),
    LegendEntry("No CBF"),
    Plot({plt_theme...,dashed}, Coordinates([0,T], [pmax, pmax])),
    [raw"\filldraw[gray, thick, opacity=0.4] (0,2) -- (10,2) -- (10,4) -- (0, 4) -- cycle;"],
    TextNode(1.3, 1.0, raw"{\footnotesize$x_{\max}$}"),
)
pgfsave("cartpole_wall_x_small.pdf", fig1)

fig2 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "thick"}, Coordinates([0, T], [θmax, θmax])),
    Plot({"smooth", "thick"}, Coordinates([0, T], -[θmax, θmax])),
)

# Plot in state space?
fig3 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$\theta$",
        xmin=-5.5,
        xmax=3.5,
        ymin=-1.0,
        ymax=1.0,
        # width="1.5in",
        # height="1.5in",
    },
    [raw"\filldraw[gray, thick, opacity=0.4] (2,1) -- (3.5,1) -- (3.5,-1) -- (2, -1) -- cycle;"],
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    Plot({plt_theme...}, Coordinates([pmax, pmax], [-2, 2])),
    # Plot({plt_theme..., dashed,}, Coordinates([-6, 4], [θd + θmax, θd + θmax])),
    # Plot({plt_theme..., dashed,}, Coordinates([-6, 4], [θd - θmax, θd - θmax])),
    TextNode(p0+0.6, θ0, raw"{\footnotesize$\mathbf{q}_0$}"),
    # TextNode(-2, -0.5, raw"{$\theta_{\max}$}"),
    TextNode(1.2, -0.7, raw"{\footnotesize$x_{\max}$}"),
    Plot({mark="*",color="black"}, Coordinates([p0], [θ0])),
)

# pgfsave("cartpole_wall.pdf", fig3)

# Plot safe set?
contour_opts = @pgf {
    contour_lua={labels=false,levels=0.0},
    line_width="1.0pt",
}

x1 = 0.0:0.1:3.0
x2 = -3.0:0.1:3.0
fig4 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$\dot{x}$",
        xmin=0.0,
        xmax=3,
        ymin=-3,
        ymax=3,
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
        Table(x1, x2, [h([x,0.0], [y,0.0]) for x in x1, y in x2]),
    ),   
    Plot({plt_theme..., dashed}, Coordinates([pmax, pmax], [-3, 3])),
    [raw"\filldraw[gray, thick, opacity=0.4] (2,-3) -- (2,3) -- (4,3) -- (4, -3) -- cycle;"],
)
# pgfsave("x_safe_set.pdf", fig4)