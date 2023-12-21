using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Colors
using Contour

# Safety constraint and gradient
h0(x) = 1 - x^2
∇h0(x) = ForwardDiff.derivative(h0, x)

# Class K function
α0(r) = r

# Top-level dynbamics
f0(x) = 0.0
g0(x) = 1.0

# Top-level CBF controller
a0(x) = α0(h0(x))
b0(x) = norm(∇h0(x))^2
λ0(a,b) = ReducedOrderModelCBFs.λSoftplus(a, b, 0.1)
k0(x) = λ0(a0(x), b0(x))*∇h0(x)

# CBF
μ = 1.0
h(x, y, μ) = h0(x) - 0.5*inv(μ)*norm(y - k0(x))^2

# Plot safe set
xs = -1.2:0.01:1.2
ys = -2.0:0.01:2.0

### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

# Plot safe set
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
    Table(contours(xs, ys, h.(xs, ys', 0.5), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[3]},
        opacity=1.0,
    },
    Table(contours(xs, ys, h.(xs, ys', 1.0), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[4]},
        opacity=1.0,
    },
    Table(contours(xs, ys, h.(xs, ys', 2.0), [0.0]))),
    Plot({
        plt_theme...,
        contour_prepared={labels=false, draw_color=colors[2]},
        opacity=0.5,
    },
    Table(contours(xs, ys, h.(xs, ys', 10.0), [0.0]))),
)

# pgfsave("double_int_safesets.pdf", ax)