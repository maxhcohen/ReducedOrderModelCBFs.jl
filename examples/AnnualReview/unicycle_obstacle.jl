using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Contour

# Get system models
Σ0 = SingleIntegrator(2); # Reduced-order model
Σ = Unicycle(); # Full-order model

# Obstacle and CBF for ROM
xo = [-1.0, 1.0]
ro = 0.4
h0(x) = norm(x - xo)^2 - ro^2
cbf0 = ControlBarrierFunction(h0, s -> s);

# Desired controller for ROM
kd0(x) = -x

# Smooth safety filter for reduced-order model
k0 = SmoothSafetyFilter(Σ0, cbf0, kd0, formula="gaussian", σ=0.1);

# Now decompose safe velocity into safe forward velocity and heading
k0v(x) = norm(k0(x))
k0ω(x) = k0(x)/norm(k0(x)) # Only valid for k0(x) ≠ 0

# Convert safe heading direction into safe heading angle
ψ0(x) = atan(k0ω(x)[2], k0ω(x)[1])

# Bottom layer state of cascaded system
ξ(x) = [cos(x[3]), sin(x[3])]

# CBF for full-order dynamics
μ = 1.0
h(x) = h0(x[1:2]) - (0.5/μ)*norm(ξ(x) - k0ω(x[1:2]))^2
cbf = ControlBarrierFunction(h, s -> s);

# Desired controller for full-order system
Kp = 0.2
Kψ = 3
qd = [0.0, 0.0]
kd(x) = [Kp*norm(x[1:2] - qd), -Kψ*(sin(x[3]) - sin(ψ0(x[1:2])))]

# Now build safety filter for full-order system
k = ReluSafetyFilter(Σ, cbf, kd);

# Initial conditions
q0 = [-2.1, 2.0]
ω0 = 0.0
x0 = [q0; ω0]

# Run simulation
T = 20.0
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