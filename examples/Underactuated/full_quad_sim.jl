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

# Don't actually need these Lie derivatives but can check to make sure relative degree holds
Lfy1(x) = ForwardDiff.gradient(y1, x)'*Σ.f(x)
Lfy2(x) = ForwardDiff.jacobian(y2, x)*Σ.f(x)
Lgy2(x) = ForwardDiff.jacobian(y2, x)*Σ.g(x)
LgLfy1(x) = ForwardDiff.gradient(Lfy1, x)'*Σ.g(x)

# Total output
y(x) = [y1(x); y2(x)]
Lfy(x) = ForwardDiff.jacobian(y, x)*Σ.f(x)

# Make CBF just based on outputs
zmin = 0.5
λ = 0.1
R33(x, y) = 1.0 - 2*x^2 - 2*y^2
# h1(y1, y2) = y1 - zmin - λ*(1.0 - ez'*rotmatrix_from_quat(quat(y2...))*ez)
h1(y1, y2) = y1 - zmin - λ*(1.0 - R33(y2...))
h1(y) = h1(y[1], y[2:3])

# Reduced-order model: 3D single integrator
Σ0 = CustomControlAffineSystem(3, 3, x -> zeros(3), x -> diagm(ones(3)))

# Smooth safety filter
α0 = 1.0
k1 = SmoothSafetyFilter(Σ0, h1, r -> α0*r, x -> zeros(3))

# Now build CBF for full-order system
μ = 1.0
h(x) = h1(y(x)) - (0.5/μ)*norm(Lfy(x)- k1(y(x)))^2

# Get nominal controller
kDF = DiffFlatQuadController(xd=[2.0, 1.0, 0.0])

# Get safety filer
α = α0
kCBF = ReluSafetyFilter(Σ, h, r -> α*r, x -> kDF(x))

# Initial states
p = [0.0, 0.0, 3.0]
q = [1.0, 0.0, 0.0, 0.0]
v = zeros(3)
x = [p; q; v]

# Simulate
T = 20.0
sol = simulate(Σ, x, kCBF, T)
ts = 0.0:0.01:T

# Plot results
ax_theme = get_ieee_theme_large()
plt_theme = get_plt_theme()
colors = get_colors()

# Positions
fig1 = @pgf Axis(
    {
        width="6in",
        height="4in",
        xlabel=raw"$t$",
        ylabel=raw"$p(t)$",
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
fig1 = @pgf Axis(
    {
        ax_theme...,
        xlabel=raw"$t$",
        xmin=ts[1],
        xmax=ts[end],
        height="1.5in",
    },
    [raw"\filldraw[gray, thick, opacity=0.4] (0,0) -- (0,-1) -- (20,-1) -- (20, 0.0) -- cycle;"],
    Plot({plt_theme..., color=colors[1],}, Coordinates(ts, h1.(y.(sol.(ts))))),
    LegendEntry(raw"$\psi(\mathbf{y}(\mathbf{x}(t)))$"),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, h.(sol.(ts)))),
    LegendEntry(raw"$h(\mathbf{x}(t))$"),
)

# pgfsave("3d_quad_barriers.pdf", fig1)