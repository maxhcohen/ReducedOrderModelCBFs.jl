using Revise
using ReducedOrderModelCBFs
using PGFPlotsX
using LaTeXStrings

# Region to plot over
a = -0.5:0.01:0.5
b = 1.0
σ = 0.1

# PGF plot it
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()
@pgf ax = Axis(
    {
        ax_theme...,
        xlabel=L"a",
        ylabel=L"\lambda(a,b)",
        xmin=-0.5,
        xmax=0.5,
        width="3.5in",
        height="2.22in",
    },
    Plot({plt_theme..., color="black"}, Coordinates(a, ReducedOrderModelCBFs.λRelu.(a, b))),
    LegendEntry("ReLU"),
    Plot({plt_theme..., color=colors[1]}, Coordinates(a, ReducedOrderModelCBFs.λSontag.(a, b, σ))),
    LegendEntry("Sontag"),
    Plot({plt_theme..., color=colors[2]}, Coordinates(a, ReducedOrderModelCBFs.λHalfSontag.(a, b, σ))),
    LegendEntry("Half-Sontag"),
    Plot({plt_theme..., color=colors[3]}, Coordinates(a, ReducedOrderModelCBFs.λSoftplus.(a, b, σ))),
    LegendEntry("Softplus"),
    Plot({plt_theme..., color=colors[4]}, Coordinates(a, ReducedOrderModelCBFs.λGaussian.(a, b, σ))),
    LegendEntry("Gaussian"),
)

# Now look at lambda for fixed a
a = 1.0
b = -2.0:0.01:2.0
σ = 0.3

@pgf ax2 = Axis(
    {
        ax_theme...,
        xlabel=L"b",
        ylabel=L"\lambda(a,b)",
        xmin=-1.0,
        xmax=2.0,
        width="3.5in",
        height="2.22in",
        legend_pos="north west",
        # title=L"a=1",
    },
    Plot({plt_theme..., color="black"}, Coordinates(b, ReducedOrderModelCBFs.λRelu.(a, b))),
    LegendEntry("ReLU"),
    Plot({plt_theme..., color=colors[1]}, Coordinates(b, ReducedOrderModelCBFs.λSontag.(a, b, σ))),
    LegendEntry("Sontag"),
    Plot({plt_theme..., color=colors[2]}, Coordinates(b, ReducedOrderModelCBFs.λHalfSontag.(a, b, σ))),
    LegendEntry("Half-Sontag"),
    Plot({plt_theme..., color=colors[3]}, Coordinates(b, ReducedOrderModelCBFs.λSoftplus.(a, b, σ))),
    LegendEntry("Softplus"),
    Plot({plt_theme..., color=colors[4]}, Coordinates(b, ReducedOrderModelCBFs.λGaussian.(a, b, σ))),
    LegendEntry("Gaussian"),
)


# Now look at lambda for fixed a
a = -1.0
b = 0.1:0.01:2.0

@pgf ax3 = Axis(
    {
        ax_theme...,
        xlabel=L"b",
        ylabel=L"\lambda(a,b)",
        xmin=0.5,
        xmax=2.0,
        width="3.33in",
        height="2.22in",
        legend_pos="north east",
        title=L"a=-1",
    },
    Plot({plt_theme..., color="black"}, Coordinates(b, ReducedOrderModelCBFs.λRelu.(a, b))),
    LegendEntry("ReLU"),
    Plot({plt_theme..., color=colors[1]}, Coordinates(b, ReducedOrderModelCBFs.λSontag.(a, b, σ))),
    LegendEntry("Sontag"),
    Plot({plt_theme..., color=colors[2]}, Coordinates(b, ReducedOrderModelCBFs.λHalfSontag.(a, b, σ))),
    LegendEntry("Half-Sontag"),
    Plot({plt_theme..., color=colors[3]}, Coordinates(b, ReducedOrderModelCBFs.λSoftplus.(a, b, σ))),
    LegendEntry("Softplus"),
    Plot({plt_theme..., color=colors[4]}, Coordinates(b, ReducedOrderModelCBFs.λGaussian.(a, b, σ))),
    LegendEntry("Gaussian"),
)

pgfsave("smooth_formulas.pdf", ax)
pgfsave("smooth_formulas_b1.pdf", ax2)
pgfsave("smooth_formulas_b3.pdf", ax3)