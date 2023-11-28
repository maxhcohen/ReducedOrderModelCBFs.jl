using Revise
using ReducedOrderModelCBFs
using LinearAlgebra

# System models
Σ0 = SingleIntegrator(2); # Reduced-order model
Σ = DoubleIntegrator(RoboticSystem, 2); # Full-order model

# Obstacle and CBF for ROM
xo = [-1.0, 1.0]
ro = 0.4
h0(x) = norm(x - xo)^2 - ro^2
cbf0 = ControlBarrierFunction(h0, s -> s);

# Desired controller for ROM
kd0(x) = -x

# Smooth safety filter for reduced-order model
k0 = SmoothSafetyFilter(Σ0, cbf0, kd0, formula="gaussian", σ=0.1);
kID = InverseDynamicsCBFQP(Σ, cbf0, k0, μ=5.0, Kp=10.0, Kd=1.0)

# Initial conditions
q0 = [-2.1, 2.0]
q̇0 = zeros(2)

# Simulate
T = 10.0
sol = simulate(Σ, kID, q0, q̇0, T)

# Plot results
using PGFPlotsX
using LaTeXStrings
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

ts = 0.0:0.01:T
@pgf ax = Axis(
    {
        ax_theme...,
        xlabel=L"x_1",
        ylabel=L"x_2",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    [raw"\filldraw[color=black, fill=black!25, thick](-1,1) circle (0.4);"],
)