using Revise
using LinearAlgebra
using ForwardDiff
using ReducedOrderModelCBFs
using PGFPlotsX

# Load in reduced-order model - single integrator with dimension 2
Σ0 = SingleIntegrator(2)

# Safety constraint for reduced-order model
xmax = 2.0
h0(x) = xmax - x[1] # Don't crash into a wall located at x = 2m
# cbf = ControlBarrierFunction(h0, s -> 0.5*s)
cbf = ControlBarrierFunction(h0, s -> 0.3*s)

# Desired input for reduced order model
kd0(x) = [1.0, 0.0] # Move forward with a desired velocity of 1 m/s

# Reduced-Order Safety filter
k0_sontag = SmoothSafetyFilter(Σ0, cbf, kd0, formula="sontag", σ=0.1);
k0_half_sontag = SmoothSafetyFilter(Σ0, cbf, kd0, formula="half-sontag", σ=0.1);
k0_softplus = SmoothSafetyFilter(Σ0, cbf, kd0, formula="softplus", σ=0.1);
k0_gaussian = SmoothSafetyFilter(Σ0, cbf, kd0, formula="gaussian", σ=0.1);

# Now load in full-order model
Σ = PlanarSegway(RoboticSystem)

# Initial states
q0 = [0.0, -0.138]
q̇0 = zeros(2)

# Create custom PD tracking controller
struct SegwayTrackingController <: FeedbackController
    get_input::Function
end
(k::SegwayTrackingController)(q,q̇) = k.get_input(q,q̇)

# Constructor
function SegwayTrackingController(Kp, Kφ, Kφ̇, k0::SmoothSafetyFilter)
    function pd_controller(q, q̇)
        p = q[1]
        φ = q[2]
        ṗ = q̇[1]
        φ̇ = q̇[2]

        return Kp*(ṗ - k0(q)[1]) + Kφ*φ + Kφ̇*φ̇
    end

    return SegwayTrackingController(pd_controller)
end

# Instantiate tracking controller
Kṗ = 50.0
# Kṗ = 30.0
Kφ = 150.0
Kφ̇ = 40.0
k_sontag = SegwayTrackingController(Kṗ, Kφ, Kφ̇, k0_sontag)
k_half_sontag = SegwayTrackingController(Kṗ, Kφ, Kφ̇, k0_half_sontag)
k_softplus = SegwayTrackingController(Kṗ, Kφ, Kφ̇, k0_softplus)
k_gaussian = SegwayTrackingController(Kṗ, Kφ, Kφ̇, k0_gaussian)

# Now sim full system
T = 10.0
ts = 0.0:0.01:T
sol_sontag = simulate(Σ, k_sontag, q0, q̇0, T)
sol_half_sontag = simulate(Σ, k_half_sontag, q0, q̇0, T)
sol_softplus = simulate(Σ, k_softplus, q0, q̇0, T)
sol_gaussian = simulate(Σ, k_gaussian, q0, q̇0, T)

### Plot setup
ax_theme = get_ax_theme()
plt_theme = get_plt_theme()
colors = get_colors()

ax = @pgf Axis(
    {
        ax_theme...,
        xmin = 0,
        xmax=T,
        ymin=-0.3,
        ymax=2.3,
        width="1.75in",
        height = "1.5in",
        tick_label_style={font="\\footnotesize"},
    },
    Plot({plt_theme..., color=colors[1]}, Coordinates(ts, sol_sontag.(ts, idxs=1))),
    Plot({plt_theme..., color=colors[2]}, Coordinates(ts, sol_half_sontag.(ts, idxs=1))),
    Plot({plt_theme..., color=colors[3]}, Coordinates(ts, sol_softplus.(ts, idxs=1))),
    Plot({plt_theme..., color=colors[4]}, Coordinates(ts, sol_gaussian.(ts, idxs=1))),
    Plot({plt_theme..., color="black", "dashed"}, Coordinates([0, T], [xmax, xmax])),
)

# pgfsave("planar_segway1.pdf", ax)
# pgfsave("planar_segway2.pdf", ax)
# pgfsave("planar_segway3.pdf", ax)
# pgfsave("planar_segway4.pdf", ax)