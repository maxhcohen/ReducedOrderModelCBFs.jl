using Revise
using LinearAlgebra
using ReducedOrderModelCBFs
using ForwardDiff
using PGFPlotsX
using LaTeXStrings
using Contour

# Get system
Σ = SingleIntegrator();

# Obstacle and CBF
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x - xo)^2 - ro^2
∇h(x) = ForwardDiff.gradient(h, x)
cbf = ControlBarrierFunction(h);

# Desired controller
kd(x) = -x

# Make cbfs and controller for differnt class K functions
cbf1 = ControlBarrierFunction(h, r -> r);
cbf2 = ControlBarrierFunction(h, r -> 0.2*r);
cbf3 = ControlBarrierFunction(h, r -> 3*r);

k1 = ReluSafetyFilter(Σ, cbf1, kd);
k2 = ReluSafetyFilter(Σ, cbf2, kd);
k3 = ReluSafetyFilter(Σ, cbf3, kd);

# simulate
x0 = [-2.1, 2.0]
T = 10.0
sol1 = simulate(Σ, k1, x0, T)
sol2 = simulate(Σ, k2, x0, T)
sol3 = simulate(Σ, k3, x0, T)
sol4 = simulate(Σ, (x, t) -> kd(x), x0, T)
ts = 0.0:0.05:T

# Plot results
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

# Plot trajectories
@pgf ax1 = Axis(
    {
        ax_theme...,
        # xlabel=L"x_1",
        # ylabel=L"x_2",
        width="3.33in",
        height="3.0in",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(sol1.(ts, idxs=1), sol1.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(sol2.(ts, idxs=1), sol2.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(sol3.(ts, idxs=1), sol3.(ts, idxs=2))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(sol4.(ts, idxs=1), sol4.(ts, idxs=2))),
    [raw"\filldraw[color=black, fill=black!25, thick](-1,1) circle (0.4);"],
)
pgfsave("single_int_states.pdf", ax1)

# Plot CBF evolution
@pgf ax2 = Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        width="3.33in",
        height="1.75in",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, h.(sol1.(ts)))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, h.(sol2.(ts)))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, h.(sol3.(ts)))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(ts, h.(sol4.(ts)))),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [0, 0])),
)
pgfsave("single_int_h.pdf", ax2)

# Plot control inputs
@pgf ax3 = Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        width="3.33in",
        height="1.75in",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, norm.(k1.(sol1.(ts))))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, norm.(k2.(sol2.(ts))))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, norm.(k3.(sol3.(ts))))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(ts, norm.(kd.(sol4.(ts))))),
)
pgfsave("single_int_controls.pdf", ax3)

# Plot rate of CBF
ḣ(x, u) = ∇h(x)'*(Σ.f(x) + Σ.g(x)*u)
@pgf ax4 = Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        width="1.75in",
        height="1.75in",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, ḣ.(sol1.(ts), k1.(sol1.(ts))))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, ḣ.(sol2.(ts), k2.(sol2.(ts))))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, ḣ.(sol3.(ts), k3.(sol3.(ts))))),
)