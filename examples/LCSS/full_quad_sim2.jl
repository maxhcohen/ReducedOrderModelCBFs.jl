using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX
using Colors

# Load in model
Σ = Quadrotor()

# Unit vector in z direction
ez = [0.0, 0.0, 1.0]

# Define some outputs
y1(x) = x[3] # z position
y2(x) = x[5:6] # Second and third components of quaternion

# Total output
y(x) = [y1(x); y2(x)]
Lfy(x) = ForwardDiff.jacobian(y, x)*Σ.f(x)

# Make CBF just based on outputs
zmin = 0.5
λ = 0.1
R33(x, y) = 1.0 - 2*x^2 - 2*y^2
ψ(y1, y2) = y1 - zmin - λ*(1.0 - R33(y2...))
ψ(y) = ψ(y[1], y[2:3])

# Reduced-order model: 3D single integrator
Σ0 = CustomControlAffineSystem(3, 3, x -> zeros(3), x -> diagm(ones(3)))

# Smooth safety filter
α0 = 1.0
k1 = SmoothSafetyFilter(Σ0, ψ, r -> α0*r, x -> zeros(3))

# Now build CBF for full-order system
μ = 1.0
h(x) = ψ(y(x)) - (0.5/μ)*norm(Lfy(x)- k1(y(x)))^2

# Get two different differential flattness controllers
xd1 = [1.0, 3.0, 4.0]
xd2 = [-2.0, 0.0, -1.0]
kDF1 = DiffFlatQuadController(xd=xd1)
kDF2 = DiffFlatQuadController(xd=xd2)

# Define time-varying nominal controller
function kd(x,t)
    if t ≥ 0 && t < 10.0
        return kDF1(x)
    elseif t ≥ 10.0 && t < 20.0
        return kDF2(x)
    elseif t ≥ 20.0 && t < 30.0
        return kDF1(x)
    elseif t ≥ 30.0 && t < 40.0
        return kDF2(x)
    elseif t ≥ 40.0 && t < 50.0
        return kDF1(x)
    elseif t ≥ 50.0 && t < 60.0
        return kDF2(x)
    elseif t ≥ 60.0
        return kDF1(x)
    end
end

# Get safety filer
α = α0
kCBF = TimeVaryingReluSafetyFilter(Σ, h, r -> α*r, kd)

# Initial states
p = [0.0, 0.0, 3.0]
q = [1.0, 0.0, 0.0, 0.0]
v = zeros(3)
x = [p; q; v]

# Simulate
T = 70.0
sol = simulate(Σ, (x,t) -> kCBF(x,t), x, T)
ts = 0.0:0.01:T

# Plot results
ax_theme = get_ieee_theme_large()
plt_theme = get_plt_theme()
colors = get_colors()

# Positions
fig1 = @pgf Axis(
    {
        ax_theme...,
        width="6in",
        height="4in",
        xlabel=raw"$t$",
        ylabel=raw"$\mathbf{p}(t)$",
        thick,
        xmin=ts[1],
        xmax=ts[end],
    },
    Plot({line_width="1.5pt", color=colorant"dodgerblue3",}, Coordinates(ts, sol.(ts, idxs=1))),
    LegendEntry(raw"$x(t)$"),
    Plot({line_width="1.5pt", color=colorant"firebrick",}, Coordinates(ts, sol.(ts, idxs=2))),
    LegendEntry(raw"$y(t)$"),
    Plot({line_width="1.5pt", color=colorant"green",}, Coordinates(ts, sol.(ts, idxs=3))),
    LegendEntry(raw"$z(t)$"),
)

# Barrier function and constraints
fig2 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$t$",
        xmin=ts[1],
        xmax=ts[end],
        height="1.5in",
    },
    [raw"\filldraw[gray, thick, opacity=0.4] (0,0) -- (0,-1) -- (70,-1) -- (70, 0.0) -- cycle;"],
    Plot({plt_theme..., color=colors[1],}, Coordinates(ts, h1.(y.(sol.(ts))))),
    LegendEntry(raw"$\psi(\mathbf{y}(\mathbf{x}(t)))$"),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, h.(sol.(ts)))),
    LegendEntry(raw"$h(\mathbf{x}(t))$"),
)

# pgfsave("3d_quad_barriers.pdf", fig2)

# Control inputs?
u(t) = kCBF(sol(t), t)
u1(t) = u(t)[1]
u2(t) = u(t)[2]
u3(t) = u(t)[3]
u4(t) = u(t)[4]

fig3 = @pgf Axis(
    {
        width="6in",
        height="4in",
        xlabel=raw"$t$",
        ylabel=raw"$u(t)$",
        thick,
        xmin=ts[1],
        xmax=ts[end],
    },
    Plot({line_width="1.5pt", color=colors[1],}, Coordinates(ts, u1.(ts))),
    LegendEntry(raw"$\tau(t)$"),
    Plot({line_width="1.5pt", color=colors[2],}, Coordinates(ts, u2.(ts))),
    LegendEntry(raw"$\omega_1(t)$"),
    Plot({line_width="1.5pt", color=colors[3],}, Coordinates(ts, u3.(ts))),
    LegendEntry(raw"$\omega_2(t)$"),
    Plot({line_width="1.5pt", color=colors[4],}, Coordinates(ts, u4.(ts))),
    LegendEntry(raw"$\omega_3(t)$"),
)

# Plot quaternion evolution
fig4 = @pgf Axis(
    {
        ax_theme...,
        width="6in",
        height="4in",
        xlabel=raw"$t$",
        ylabel=raw"$\mathbf{q}(t)$",
        thick,
        xmin=ts[1],
        xmax=ts[end],
    },
    Plot({line_width="1.5pt", color=colors[1],}, Coordinates(ts,  sol.(ts, idxs=4))),
    LegendEntry(raw"$q_w(t)$"),
    Plot({line_width="1.5pt", color=colors[2],}, Coordinates(ts,  sol.(ts, idxs=5))),
    LegendEntry(raw"$q_x(t)$"),
    Plot({line_width="1.5pt", color=colors[3],}, Coordinates(ts,  sol.(ts, idxs=6))),
    LegendEntry(raw"$q_y(t)$"),
    Plot({line_width="1.5pt", color=colors[4],}, Coordinates(ts,  sol.(ts, idxs=7))),
    LegendEntry(raw"$q_z(t)$"),
)