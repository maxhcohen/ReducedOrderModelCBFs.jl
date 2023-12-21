using LinearAlgebra
using Revise
using ReducedOrderModelCBFs
using PGFPlotsX
using ForwardDiff

# Make Pendubot
Σ = CartPole(RoboticSystem)

# Get inertia matrix
D(q) = Σ.M(q)

# Get actuation matrix
B = Σ.B(rand(2))

# Left anhilator of matrix
Bp = [0.0, 1.0]'

# Decouping
A(q) = Bp*inv(D(q))*B
A(q1, q2) = A([q1, q2])

# Now make CBF for system
θmax = π/6
h0(q) = θmax^2 - (q[2] - π)^2
∇h0(q) = ForwardDiff.gradient(h0, q)
a0(q) = h0(q)
b0(q) = norm(∇h0(q))^2
k0(q) = λSoftplus(a0(q), b0(q), 0.1)*∇h0(q)
h(q, q̇) = h0(q) - 0.5*(q̇[2] - k0(q)[2])^2
h(x) = h(x[1:2], x[3:4])
∇h(x) = ForwardDiff.gradient(h, x)

# Now get control affine representation
Σ2 = CartPole(ControlAffineSystem)
f(x) = Σ2.f(x)
g(x) = Σ2.g(x)
Lfh(x) = ∇h(x)'f(x)
Lgh(x) = dot(∇h(x), g(x))

# Desired controller - track a angle reference
function θd(t)
    if t < 5
        return π - π/4
    else
        return π + π/4
    end
end
kd(x, t) = -50*(x[2] - θd(t)) - 30*x[4]

# Safety filter
a(x, t) = Lfh(x) + Lgh(x)*kd(x, t) + h(x)
b(x) = norm(Lgh(x))^2
kqp(x, t) = kd(x, t) + λRelu(a(x, t), b(x))*Lgh(x)'

# Initial conditions
p0 = 0.0
θ0 = π
ṗ0 = 0.0
θ̇0 = 0.0
q0 = [p0, θ0]
q̇0 = [ṗ0, θ̇0]
x0 = vcat(q0, q̇0)

# Simulation for sanity check
T = 10.0
ts = 0.0:0.01:T
sol = simulate(Σ2, kqp, x0, T)

@pgf Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        height="1.75in",
        width="1.75in",
        # xlabel=L"t",
        # ylabel=L"\theta(t)",
    }, 
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [π+θmax, π+θmax])),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [π-θmax, π-θmax])),
    Plot({plt_theme..., color=colors[2], opacity=0.5}, Coordinates([0, 5], [π-π/4, π-π/4])),
    Plot({plt_theme..., color=colors[2], opacity=0.5}, Coordinates([5, T], [π+π/4, π+π/4])),
)

@pgf Axis(
    {
        ax_theme...,
        xmin=0,
        xmax=T,
        width="1.75in",
        height="1.75in",
    }, 
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, h0.(sol.(ts, idxs=1:2)))),
    Plot({plt_theme..., "black", "dashed"}, Coordinates([0,T], [0, 0])),
)

# Plot surface of A
xs = range(-π/2, π/2, 100)
ys = range(-π/2, π/2, 100)
### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

@pgf Axis(
    {
        ax_theme...,
        colorbar,
        "colormap/viridis",
    },
    Plot3(
        {
            surf,
            shader = "flat",
        },
        Table(xs, ys, [A([x, y]) for x in xs, y in ys])
    )
)

@pgf Axis(
    {
        ax_theme...,
        xmin=-π/2,
        xmax=π/2,
        xlabel=raw"$q_1$",
        ylabel=raw"$A(q)$",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(xs, A.(xs, 0.0))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(xs, A.(xs, π/4))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(xs, A.(xs, π/2))),
    Plot({plt_theme...,}, Coordinates([-π/2, π/2], [0.0, 0.0])),
)

@pgf Axis(
    {
        ax_theme...,
        xmin=-π/2,
        xmax=π/2,
        xlabel=raw"$q_2$",
        ylabel=raw"$A(q)$",
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ys, A.(0.0, ys))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ys, A.(π/4, ys))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(ys, A.(π/2, ys))),
    Plot({plt_theme...,}, Coordinates([-π/2, π/2], [0.0, 0.0])),
)