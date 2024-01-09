using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in model
Σr = CartPole(RoboticSystem)

# System output - pendulum angle
y(q) = q[2]
J(q) = [0.0, 1.0]'

# Desired pendulum angle and maximum deviation
θd = 0.0
θmax = π/12

# Output constraint
h0(y) = θmax^2 - (θd - y)^2

# Desired controller for full-order system: move to a desired location
xd = 3.0
Kp = 1.0
Kd = 2.0
kd1(q, q̇) = -Kp*(q[1] - xd) - Kd*q̇[1] + 0.05*q̇[2]
# kd1(q, q̇) = 0.0

# Desired controller for reduced-order system: keep angle at θd
kd0(y) = -y

# Reduced-order CBF parameters
α = 5.0
∇h0(y) = ForwardDiff.derivative(h0, y)
Lfh0(y) = 0.0
Lgh0(y) = ∇h0(y)
a(y) = Lfh0(y) + Lgh0(y)*kd0(y) + α*h0(y)
b(y) = norm(Lgh0(y))^2

# Smooth safety filter
σ = 0.05
k0(q) = λSoftplus(a(y(q)), b(y(q)), σ)*Lgh0(y(q))'

# Now build CBF for full-order system
μ = 10.0
h(q,q̇) = h0(y(q)) - (0.5/μ)*norm(J(q)*q̇ - k0(q))^2

# Get Lie derivatives and stuff
n, m, f, g = ReducedOrderModelCBFs.to_control_affine(Σr)
dhdq(q, q̇) = ForwardDiff.gradient(q -> h(q,q̇), q)'
dhdq̇(q, q̇) = ForwardDiff.gradient(q̇ -> h(q,q̇), q̇)'
Dh(q, q̇) = hcat(dhdq(q, q̇), dhdq̇(q, q̇)) 
Lfh(q, q̇) = Dh(q, q̇)*f(vcat(q,q̇))
Lgh(q, q̇) = Dh(q, q̇)*g(vcat(q,q̇))

# Safety filter
k(q, q̇) = kd1(q,q̇) + λRelu(Lfh(q,q̇) + Lgh(q,q̇)*kd1(q,q̇) + α*h(q,q̇), norm(Lgh(q,q̇))^2)*Lgh(q,q̇)'

# Check out decoupling matrix
D = Σr.M
B = Σr.B
A(q) = J(q)*inv(D(q))*B(q)

# Initial conditions
p0 = -5.0
θ0 = 0.2
ṗ0 = 0.0
θ̇0 = 0.0
q0 = [p0, θ0]
q̇0 = [ṗ0, θ̇0]

# Sim parameters
T = 10.0
dt = 0.05
ts = 0.0:dt:T

# Simulate
sol = simulate(Σr, (q,q̇,t) -> k(q,q̇), q0, q̇0, T)

# Plot results
fig1 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$p(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=1))),
)

fig2 = @pgf Axis(
    {
        xlabel=raw"$t$",
        ylabel=raw"$\theta(t)$",
    },
    Plot({"smooth", "thick"}, Coordinates(ts, sol.(ts, idxs=2))),
    Plot({"smooth", "thick"}, Coordinates([0, T], [θmax, θmax])),
    Plot({"smooth", "thick"}, Coordinates([0, T], -[θmax, θmax])),
)

fig3 = @pgf Axis(
    {
        xlabel=raw"$p$",
        ylabel=raw"$\theta$",
    },
    Plot({"smooth", "thick", color="blue"}, Coordinates(sol.(ts, idxs=1), sol.(ts, idxs=2))),
    Plot({"smooth", "thick"}, Coordinates([pmax, pmax], [-θmax, θmax])),
    Plot({"smooth", "thick"}, Coordinates([-5, pmax], [θmax, θmax])),
    Plot({"smooth", "thick"}, Coordinates([-5, pmax], -[θmax, θmax])),
)