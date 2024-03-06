using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX
using Colors

# Dynamics
f(x) = [x[2], 0.0]
g(x) = [0.0, 1.0]

# Output
y(x) = x[1]

# Output dynamics
Lfy(x) = x[2]

# Constraint function
ymin = 0.0
ψ(y) = y - ymin

# Reduced-order model: 1D single integrator
Σ0 = CustomControlAffineSystem(1, 1, x -> 0.0, x -> 1.0)

# Smooth safety filter
α0 = 1.0
k = SmoothSafetyFilter(Σ0, ψ, r -> α0*r, x -> 0.0)

# CBF
h(x) = ψ(y(x)) - 0.5*(Lfy(x) - k(y(x)))^2

# Something that looks like the CBF
μ = 1.5
cbf1(x) = x == 0.0 ? 0.0 : μ*sqrt(x)
cbf2(x) = x == 0.0 ? 0.0 : -μ*sqrt(x)
cbf_x = 0.0:0.01:2


# Region to plot over
x1 = -2.0:0.05:2.0
x2 = -2.0:0.05:2.0


# Make plot
colors = get_colors()
contour_opts = @pgf {
    contour_lua={labels=false,levels=0.0},
    line_width="1.0pt",
}
ax_theme = get_ieee_theme_small()
fig1 = @pgf Axis(
    {
        ax_theme...,
        width="1.5in",
        height="1.5in",
        xmin=-1,
        xmax=1,
        ymin=-2,
        ymax=2,
        xticklabel="\\empty",
        yticklabel="\\empty",
        # major_tick_length="0.075cm",
        view="{0}{90}",
        "colormap/greenyellow",
    },
    Plot({smooth, line_width="1.0pt",}, Coordinates([-0.0, -0.0], [-2,2])),
    Plot({smooth, line_width="1.0pt", color=colors[1],}, Coordinates(cbf_x, cbf1.(cbf_x))),
    Plot({smooth, line_width="1.0pt", color=colors[1],}, Coordinates(cbf_x, cbf2.(cbf_x))),
    # [raw"\filldraw[black, thick, opacity=0.4] (0,-2) -- (-3,-2) -- (-3,2) -- (0, 2) -- cycle;"],
    # Plot3(
    #     {
    #         contour_opts...,
    #     },
    #     Table(x1, x2, [h([x,y]) for x in x1, y in x2]),
    # )
)

# pgfsave("constraint_set.pdf", fig1)