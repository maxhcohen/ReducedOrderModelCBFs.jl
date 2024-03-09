using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σ = PlanarQuadrotor(RoboticSystem)

y(q) = [q[2], q[3]]
J(q) = ForwardDiff.jacobian(y, q)

# Constraint on output
zc = 2.0
zmax = 0.1
θmax = π/4
ellipse(x, y, xo, yo, a, b) = 1 - ((x - xo)^2/a^2) - ((y- yo)^2/b^2)
h0(y) = ellipse(y[1], y[2], zc, 0.0, zc-zmax, θmax)

# Reduced-order model
Σ0 = CustomControlAffineSystem(2, 2, x -> [0.0, 0.0], x -> diagm(ones(2)))

# Smooth safety filter
α0 = 0.1
k0 = SmoothSafetyFilter(Σ0, h0, r -> α0*r, x -> zeros(2))

# Now build CBF for full-order system
μ = 10.0
h(q,q̇) = h0(y(q)) - (0.5/μ)*norm(J(q)*q̇ - k0(y(q)))^2

# Nominal stabilizing controller
k_theta(q, q̇) = [-5*q[1] - 10*q̇[1], q[3] + 1*q̇[3]]
k_x(q, q̇) = [-q[1] - 2*q̇[1], 5*q[1] + 10*q̇[1]]

# Get safety filer
α = 1.0
kCBF = ReluSafetyFilter(Σ, h, r -> α*r, (q,q̇) -> k_x(q, q̇))

# Simulate from initial condition
x0 = 0.0
z0 = zc
θ0 = 0.5
q0 = [x0, z0, θ0]
q̇0 = zeros(3)
T = 10.0
sol = simulate(Σ, kCBF, q0, q̇0, T)


# Do some plots
ax_theme = get_ax_theme()
colors = get_colors()
ts = 0.0:0.05:T
fig1 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$t$",
        ylabel=raw"$z(t)$",
        width="3.33in",
        height="2.9in",
        xmin=0,
        xmax=T,
    },
    Plot({"smooth", very_thick}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "dashed"}, Coordinates([0, T], [0.1, 0.1]))
)

# pgfsave("quadrotor_falling.pdf", fig1)

fig2 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
        width="6in",
        height="4in",
        xmin=0,
        xmax=T,
    },
    Plot({"smooth", very_thick}, Coordinates(ts, sol.(ts, idxs=3)))
)

# Plot ellipse
ellipse(x, y, xo, yo, a, b) = 1 - ((x - xo)^2/a^2) - ((y- yo)^2/b^2)
xs = 0.0:0.1:4.0
ys = range(-π/4, π/4, 100)

fig3 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$z$",
        ylabel=raw"$\theta$",
        width="3.33in",
        height="2.9in",
        view="{0}{90}",
    },
    Plot3(
        {
            contour_lua={labels=false, levels=0.0},
        },
        Table(xs, ys, [ellipse(x, y, z0, 0.0, z0-zmax, π/6) for x in xs, y in ys])
    )
    # Plot({"smooth", very_thick}, Coordinates(sol.(ts, idxs=2), sol.(ts, idxs=3))),
    # Plot({"smooth", "dashed"}, Coordinates([0, T], [0.1, 0.1]))
)

# Plot planar position of quadrotor
l = 0.3
function plot_quad(t)
    p1 = @pgf Plot(
        {
            "smooth", 
            very_thick,
            color="black",
        },
        Coordinates(
            [
                sol.(t, idxs=1), 
                sol.(t, idxs=1) + l*cos(sol.(t, idxs=3)), 
                sol.(t, idxs=1) - l*cos(sol.(t, idxs=3)),
            ], 
            [
                sol.(t, idxs=2), 
                sol.(t, idxs=2) + l*sin(sol.(t, idxs=3)), 
                sol.(t, idxs=2) - l*sin(sol.(t, idxs=3)),
            ]
        ),
    )

    p2 = @pgf Plot({mark="*", color="black", mark_size=1.5}, Coordinates([sol.(t, idxs=1)], [sol.(t, idxs=2)]))

    return [p1, p2]
end

ax_theme = get_ieee_theme_large()
fig4 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$z$",
        xmin=-3.3,
        xmax=5,
        ymax=2.5,
        ymin=-0.5,
        height="1.8in",
    },
    Plot({"smooth", very_thick, color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    Plot({"smooth", thick, "dashed"}, Coordinates([-5, 5], [zmax, zmax])),
    plot_quad(0.0),
    plot_quad(0.3),
    plot_quad(0.5),
    plot_quad(0.75),
    plot_quad(1.0),
    plot_quad(2.0),
    plot_quad(4.0),
    plot_quad(5.0),
    plot_quad(8.0),
    plot_quad(T),
    [raw"\filldraw[gray, thick, opacity=0.4] (-4,-2) -- (5,-2) -- (5,0.1) -- (-4, 0.1) -- cycle;"],
    [raw"\node at ", Coordinate(x0-0.2, z0+0.2),raw"{\footnotesize$\mathbf{q}_0$};"],
    [raw"\node at ", Coordinate(sol(T, idxs=1), sol(T, idxs=2)-0.2),raw"{\footnotesize$\mathbf{q}_{f}$};"],
    [raw"\node at ", Coordinate(-2.95, zmax+0.15),raw"{\footnotesize$z_{\min}$};"],
    # [raw"\node at ", Coordinate(2.0, -0.2),raw"{$h_1(\mathbf{y}(\mathbf{q}))<0$};"],
    # TextNode(3.0, 1.0, raw"{$h_1(\mathbf{y}(\mathbf{q}))\geq 0$}"),
)

# pgfsave("quadrotor_falling.pdf", fig4)

# Make group plot of height and theta
ax_theme_small = get_ieee_theme_small()

fig5 = @pgf Axis(
    {
        ax_theme_small...,
        # xlabel=raw"$t$",
        ylabel=raw"$z(t)$",
        xmin=0.0,
        xmax=T,
        ymin=-0.5,
        xticklabel=raw"\empty",
    },
    Plot({"smooth", very_thick, color=colors[3]}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "dashed"}, Coordinates([0, T], [0.1, 0.1])),
)

fig6 = @pgf Axis(
    {
        ax_theme_small...,
        # xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
        xmin=0.0,
        xmax=T,
        ymin=-π/3,
        ymax=π/3,
        xticklabel=raw"\empty",
        ytick=[-1.0, 0.0, 1.0]
    },
    Plot({"smooth", very_thick, color=colors[4]}, Coordinates(ts, sol.(ts, idxs=3))),
    Plot({"smooth", "dashed"}, Coordinates([0, T], [θmax, θmax])),
    Plot({"smooth", "dashed"}, Coordinates([0, T], [-θmax, -θmax])),
)

fig7 = @pgf Axis(
    {
        ax_theme_small...,
        xlabel=raw"$t$",
        # ylabel=raw"$h(t)$",
        xmin=0.0,
        xmax=T,
    },
    Plot({"smooth", very_thick, color=colors[1],}, Coordinates(ts, h0.(y.(sol.(ts, idxs=1:3))))),
    LegendEntry(raw"$\psi$"),
    Plot({"smooth", very_thick, color=colors[2], dashed}, Coordinates(ts, h.(sol.(ts, idxs=1:3),sol.(ts, idxs=4:6)))),
    LegendEntry(raw"$h$"),
)

u(t) = kCBF(sol(t, idxs=1:3), sol(t, idxs=4:6))
u1(t) = u(t)[1]
u2(t) = u(t)[2]
fig8 = @pgf Axis(
    {
        ax_theme_small...,
        xlabel=raw"$t$",
        ylabel=raw"$\mathbf{u}(t)$",
        xmin=0.0,
        xmax=T,
    },
    Plot({"smooth", very_thick, color=colors[1]}, Coordinates(ts, u1.(ts))),
    LegendEntry(raw"$F$"),
    Plot({"smooth", very_thick, color=colors[2]}, Coordinates(ts, u2.(ts))),
    LegendEntry(raw"$M$"),
)

# Make group plot
@pgf gp = GroupPlot(
    {
        group_style = { group_size = "2 by 2", horizontal_sep="0.5in", vertical_sep="0.1in"},
    }
)

@pgf push!(gp, fig5)
@pgf push!(gp, fig6)
@pgf push!(gp, fig7)
@pgf push!(gp, fig8)

# pgfsave("quadrotor_z_theta.pdf", gp)
