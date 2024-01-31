using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σ = PlanarQuadrotor(RoboticSystem)

# Output
# y(q) = q[2]
# J(q) = ForwardDiff.gradient(y, q)'

y(q) = [q[2], q[3]]
J(q) = ForwardDiff.jacobian(y, q)

# Constraint on output
zc = 2.0
zmax = 0.1
θmax = π/4
ellipse(x, y, xo, yo, a, b) = 1 - ((x - xo)^2/a^2) - ((y- yo)^2/b^2)
h0(y) = ellipse(y[1], y[2], zc, 0.0, zc-zmax, θmax)
# h0(y) = y - zmax

# Reduced-order model: 1D single integrator
# Σ0 = CustomControlAffineSystem(1, 1, x -> 0.0, x -> 1.0)
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
l = 0.5
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

    p2 = @pgf Plot({mark="*", color="black"}, Coordinates([sol.(t, idxs=1)], [sol.(t, idxs=2)]))

    return [p1, p2]
end
fig4 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$x$",
        ylabel=raw"$z$",
        width="3.5in",
        height="3.0in",
        xmin=-5,
        xmax=5,
    },
    Plot({"smooth", very_thick}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
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
)
# pgfsave("quadrotor_falling.pdf", fig4)
# fig5 = @pgf Axis(
#     {
#         xlabel=raw"$t$",
#         ylabel=raw"$x(t)$",
#         width="6in",
#         height="4in",
#         xmin=0,
#         xmax=T,
#     },
#     Plot({"smooth", very_thick}, Coordinates(ts, sol.(ts, idxs=1)))
# )
