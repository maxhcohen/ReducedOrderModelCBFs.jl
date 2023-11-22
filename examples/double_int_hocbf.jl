using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Colors
using Contour

# Dynamics
f(x) = [x[2], 0.0]
g(x) = [0.0, 1.0]

# Safety constriant
h(x) = 1 - x[1]^2
∇h(x) = ForwardDiff.gradient(h, x)
ḣ(x) = ∇h(x)'f(x)

# Extended CBF
p1 = 0.5
p2 = 1.0
p3 = 2.0
p4 = 5.0
h1(x) = ḣ(x) + p1*h(x)
h2(x) = ḣ(x) + p2*h(x)
h3(x) = ḣ(x) + p3*h(x)
h4(x) = ḣ(x) + p4*h(x)
h1(x1, x2) = h1([x1,x2])
h2(x1, x2) = h2([x1,x2])
h3(x1, x2) = h3([x1,x2])
h4(x1, x2) = h4([x1,x2])

# Overall safe set
H1(x) = min(h(x), h1(x))
H2(x) = min(h(x), h2(x))
H3(x) = min(h(x), h3(x))
H4(x) = min(h(x), h4(x))
H1(x1, x2) = H1([x1, x2])
H2(x1, x2) = H2([x1, x2])
H3(x1, x2) = H3([x1, x2])
H4(x1, x2) = H4([x1, x2])

# Plot safe set
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()
xs = -1.2:0.01:1.2
ys = -4.0:0.01:4.0
@pgf ax = Axis(
    {
        ax_theme...,
        xlabel=L"x",
        ylabel=L"v",
        xmin=-1.2,
        xmax=1.2,
        ymin=-2.0,
        ymax=2.0,
        width="3.33in",
        height="2.22in",
    },
    Plot({plt_theme...,}, Coordinates([1.0, 1.0], [-2,2])),
    Plot({plt_theme...,}, Coordinates([-1.0, -1.0], [-2,2])),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[1]},
        opacity=1.0,
    },
    Table(contours(xs, ys, H1.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[2]},
        opacity=1.0,
    },
    Table(contours(xs, ys, H2.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[3]},
        opacity=1.0,
    },
    Table(contours(xs, ys, H3.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[4]},
        opacity=1.0,
    },
    Table(contours(xs, ys, H4.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[1]},
        opacity=0.2,
    },
    Table(contours(xs, ys, h1.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[2]},
        opacity=0.2,
    },
    Table(contours(xs, ys, h2.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[3]},
        opacity=0.2,
    },
    Table(contours(xs, ys, h3.(xs, ys'), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[4]},
        opacity=0.2,
    },
    Table(contours(xs, ys, h4.(xs, ys'), [0.0]))),
)

# pgfsave("double_int_hocbf.pdf", ax)