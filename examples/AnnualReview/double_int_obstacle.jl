using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Contour

# Get system models
Σ0 = SingleIntegrator(2); # Reduced-order model
Σ = DoubleIntegrator(ControlAffineSystem, 2); # Full-order model

# Obstacle and CBF for ROM
xo = [-1.0, 1.0]
ro = 0.4
h0(x) = norm(x - xo)^2 - ro^2
cbf0 = ControlBarrierFunction(h0, s -> s);

# Desired controller for ROM
kd0(x) = -x

# Smooth safety filter for reduced-order model
k0 = SmoothSafetyFilter(Σ0, cbf0, kd0, formula="gaussian", σ=0.1);

# CBF for full-order dynamics
μ = 5.0
h(x) = h0(x[1:2]) - (0.5/μ)*norm(x[3:4] - k0(x[1:2]))^2
cbf = ControlBarrierFunction(h, s -> s);

# Use CBF to construct QP-based safety filter
kd(x) = -x[1:2] - 2*x[3:4]
k = ReluSafetyFilter(Σ, cbf, kd);

# Initial conditions
x0 = [-2.1, 2.0, 0.0, 0.0]

# Run simulation
T = 15.0
dt = 0.01
ts = 0.0:dt:T
sol = simulate(Σ, k, x0, T)

# Plot results
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

@pgf ax = Axis(
    {
        ax_theme...,
        xlabel=L"x_1",
        ylabel=L"x_2",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    [raw"\filldraw[color=black, fill=black!25, thick](-1,1) circle (0.4);"],
)