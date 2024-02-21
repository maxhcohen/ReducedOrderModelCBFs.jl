using LinearAlgebra
using ReducedOrderModelCBFs
using PGFPlotsX
using Colors
using LaTeXStrings

# Get system
Σ = SingleIntegrator()

# Obstacle and CBF
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x - xo)^2 - ro^2
cbf = ControlBarrierFunction(h);

# Make controllers
kd(x) = -x
kQP = ReluSafetyFilter(Σ, cbf, kd);
kSontag = SmoothSafetyFilter(Σ, cbf, kd, formula="sontag", σ=0.1);
kHalfSontag = SmoothSafetyFilter(Σ, cbf, kd, formula="half-sontag", σ=0.1);
kSoftplus = SmoothSafetyFilter(Σ, cbf, kd, formula="softplus", σ=0.1);
kGaussian = SmoothSafetyFilter(Σ, cbf, kd, formula="gaussian", σ=0.1);

# Initial condition
x0 = [-2.1, 2.0]
T = 15.0

# Run sims
solQP = simulate(Σ, kQP, x0, T)
solSontag = simulate(Σ, kSontag, x0, T) # Blows up
solHalfSontag = simulate(Σ, kHalfSontag, x0, T)
solSoftplus = simulate(Σ, kSoftplus, x0, T)
solGaussian = simulate(Σ, kGaussian, x0, T)

### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

### Make plot
ts = 0.0:0.01:T
@pgf ax = Axis(
    {
        ax_theme...,
        xlabel=L"x_1",
        ylabel=L"x_2",
        xmin=-2.2,
        xmax=0.2,
        ymin=-0.2,
        ymax=2.2,
    },
    # Add a circle plot
    [raw"\filldraw[color=black, fill=black!25, thick](-1,1) circle (0.4);"],
    Plot({plt_theme..., color="black"}, Coordinates(solQP.(ts, idxs=1), solQP.(ts, idxs=2))),
    # Plot({plt_theme..., color=colors[1]}, Coordinates(solSontag.(ts, idxs=1), solSontag.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(solHalfSontag.(ts, idxs=1), solHalfSontag.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(solSoftplus.(ts, idxs=1), solSoftplus.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(solGaussian.(ts, idxs=1), solGaussian.(ts, idxs=2))),
)