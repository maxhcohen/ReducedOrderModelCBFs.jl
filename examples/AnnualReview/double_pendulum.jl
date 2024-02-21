using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX
using LaTeXStrings
using Colors

# Load in system
l1 = 1.0
l2 = 1.0
Σ = DoublePendulum(RoboticSystem, l1=l1, l2=l2)

# Function to get position of end-effector
p1(q) = l1*[sin(q[1]), -cos(q[1])]
p2(q) = p1(q) + l2*[sin(q[1] + q[2]), -cos(q[1] + q[2])]

# Want position of the end-effector within some bound
pmax = 1.5
h0(q) = pmax^2 - (p2(q)[1])^2

# Use this function as a CBF for a ROM - single integrator
kd0(q) = zeros(2)
∇h0(q) = ForwardDiff.gradient(h0, q)
a0(q) = ∇h0(q)'kd0(q) + h0(q)
b0(q) = norm(∇h0(q))^2
k0(q) = kd0(q) + λSoftplus(a0(q), b0(q), 0.1)*∇h0(q)

# Now extend this to a CBF for the overall system
μ = 10.0
h(q, q̇) = h0(q) - (0.5/μ)*(q̇ - k0(q))'*Σ.M(q)*(q̇ - k0(q))

# Now convert our guy into a control affine system
Σ1 = DoublePendulum(ControlAffineSystem)

# Make CBF
kd(x) = -1.2*x[3:4]
h(x) = h(x[1:2], x[3:4])
cbf = ControlBarrierFunction(h, s -> 5*s)
k = ReluSafetyFilter(Σ1, cbf, kd)

# Now simulate system
q0 = [π+0.1, 0.0]
q̇0 = zeros(2)
x0 = vcat(q0, q̇0)
T = 15.0
sol = simulate(Σ1, k, x0, T)
q(t) = sol(t, idxs=1:2)

### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

plot_link(t) = [[0.0, p1(q(t))[1], p2(q(t))[1]], [0.0, p1(q(t))[2], p2(q(t))[2]]]

# Make plot
@pgf ax = Axis(
    {
        ax_theme...,
        # xlabel=L"q_1",
        # ylabel=L"q_2",
        width="1.75in",
        height="1.5in",
        ymin=-2.4,
        ymax=2.4,
        xmin=-1.9,
        xmax=1.9,
    },
    Plot({plt_theme..., color=colors[2]}, Coordinates([-1.5, -1.5], [-2.5, 2.5])),
    Plot({plt_theme..., color=colors[2]}, Coordinates([1.5, 1.5], [-2.5, 2.5])),
)

ts = 0.0:0.5:6
opacities = range(1.0, 0.1, length(ts))
for (i,t) in enumerate(ts)
    @pgf push!(ax, Plot({plt_theme..., color="gray", opacity=opacities[i]}, Coordinates(plot_link(t)[1], plot_link(t)[2])))
end

# pgfsave("double_pendulum.pdf", ax)

# Make plot for safety constraint
@pgf ax2 = Axis(
    {
        ax_theme...,
        # xlabel=L"q_1",
        # ylabel=L"q_2",
        width="1.75in",
        height="1.5in",
        xmin=0.0,
        xmax=T,
        ymin=-0.5,
    },
    Plot({plt_theme..., color=colors[4]}, Coordinates(0.0:0.01:T, h0.(q.(0.0:0.01:T)))),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0.0, T], [0.0, 0.0])),
)

# pgfsave("double_pendulum_h.pdf", ax2)