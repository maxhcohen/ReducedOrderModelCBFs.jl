using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Colors
using Contour

# Open or closed-loop sim
closed_loop = true

# CBF
θ̄ = π/4
h0(x) = θ̄^2 - x^2
h(x, y) = h0(x) - 0.5*(y + x)^2
h(x) = h(x[1], x[2])
cbf = ControlBarrierFunction(h, r -> r)

# Dynamics
Σ = InvertedPendulum(ControlAffineSystem, mass=2.0, length=1.0)

# Controller
kQP = ReluSafetyFilter(Σ, cbf)

# Closed-loop vector field for plotting - choose closed
if closed_loop
    fcl(x) = Σ.f(x) + Σ.g(x)*kQP(x)
else
    fcl(x) = Σ.f(x)
end
# Grid for plotting safe set
xs = -2.0:0.01:2.0
ys = -2.0:0.01:2.0

# Grid for plotting vector field
Xs = -2.0:0.2:2.0
Ys = -2.0:0.2:2.0
XX, YY = meshgrid(Xs, Ys)

# Compute vector field over meshgrid
scale = 0.125 #0.075
normalize_arrows = true
dx, dy = mesh_vector_field(XX, YY, fcl, scale, normalize_arrows)

# Map magnitude to colors
c = norm.(fcl.(XX, YY))

### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

@pgf ax = Axis(
    {
        ax_theme...,
        xmin = -1.2,
        xmax = 1.2,
        ymin=-2,
        ymax=2,
        # xlabel=L"\theta",
        # ylabel=L"\dot{\theta}",
        # title= closed_loop ? "Closed-loop" : "Open-loop",
        "colormap/viridis",
        width="2.5in",
        height="2.0in",
    },
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[2]},
        opacity=1.0,
    },
    Table(contours(xs, ys, h.(xs, ys'), [0.0]))),
    Plot({plt_theme...,}, Coordinates([π/4, π/4], [-2,2])),
    Plot({plt_theme...,}, Coordinates([-π/4, -π/4], [-2,2])),
    Plot(
        {
            quiver = {u = "\\thisrow{u}", v = "\\thisrow{v}"},
            "-stealth",
            point_meta="\\thisrow{C}",
            "quiver/colored",
        },
        Table({meta="C"}, x=XX, y=YY, u=dx, v=dy, C=c)),
)

pgfsave("inv_pen_closed.pdf", ax)