using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX
using Colors

# Load in model
Σ = Quadrotor()

# Rotation matrix from quaternion
ez = [0.0, 0.0, 1.0]

# Define some outputs
y1(x) = x[3]
y2(x) = x[4:7]
Lfy1(x) = ForwardDiff.gradient(y1, x)'*Σ.f(x)
Lfy2(x) = ForwardDiff.jacobian(y2, x)*Σ.f(x)
Lgy2(x) = ForwardDiff.jacobian(y2, x)*Σ.g(x)

LgLfy1(x) = ForwardDiff.gradient(Lfy1, x)'*Σ.g(x)

# Make CBF just based on y1
zmin = 0.5
h1(y) = y - zmin

# Reduced-order model: 1D single integrator
Σ0 = CustomControlAffineSystem(1, 1, x -> 0.0, x -> 1.0)

# Smooth safety filter
α0 = 1.0
k1 = SmoothSafetyFilter(Σ0, h1, r -> α0*r, x ->0.0)

# Now build CBF for full-order system
μ = 1.0

# CBF with orientation penalty
h(x) = h1(y1(x)) - (0.5/μ)*norm(Lfy1(x)- k1(y1(x)))^2 - λ*(1.0 - ez'*rotmatrix_from_quat(quat(y2(x)...))*ez)

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
        width="6in",
        height="4in",
        xlabel=raw"$t$",
        ylabel=raw"Barrier value",
        thick,
        xmin=ts[1],
        xmax=ts[end],
    },
    Plot({line_width="1.5pt", color=colorant"dodgerblue3",}, Coordinates(ts, h1.(y1.(sol.(ts))))),
    LegendEntry(raw"$h_1(x(t))$"),
    Plot({line_width="1.5pt", color=colorant"purple3",}, Coordinates(ts, h.(sol.(ts)))),
    LegendEntry(raw"$h(x(t))$"),
)