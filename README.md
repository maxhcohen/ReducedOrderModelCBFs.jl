# ReducedOrderModelCBFs

This repository contains a Julia package that implements various control barrier function (CBFs) techniques based on reduced-order models. The code here is based on a collection of papers (see related works at the end of this readme) and may be used to reconstruct many of the examples from our tutorial/survey paper:

 M. H. Cohen, T. G. Molnar, and A. D. Ames, "Safety-Critical Control of Autonomous Systems: Control Barrier Functions via Reduced-Order Models," under review.

If you find this code useful, please consider citing some of the works mentioned at the end of this README.

## Installation
The most reliable way to use this code would be to clone this repo, activate the corresponding Julia environment, and then run the code examples. Alternatively, you can directly add this repo as an unregistered Julia package by entering the Julia package manager and executing:

    add https://github.com/maxhcohen/ReducedOrderModelCBFs.jl

## Quick tutorial
We'll walk step-by-step through one of the examples that demonstrates the procedure used to construct CBFs via reduced-order models. The objective here is to design a controller for a double integrator to avoid an obstacle by extending a CBF for a single integrator to that for a double integrator

First, we need to load the packages we'll use:

```julia
using LinearAlgebra
using ReducedOrderModelCBFs
```

Next, we create our reduced-order model -- a 2D single integrator:
```julia
Σ0 = SingleIntegrator(2)
```

Now, we construct a CBF for this single integrator:
```julia
xo = [-1.0, 1.0] # Obstacle center
ro = 0.4 # Obstacle radius
h0(x) = norm(x - xo)^2 - ro^2 # Distance to obstacle
cbf0 = ControlBarrierFunction(h0, s -> s); # CBF for reduced-order model
```

We use this CBF to construct a smooth safety filter for the reduced-order model:
```julia
kd0(x) = -x # Desired controller for ROM
k0 = SmoothSafetyFilter(Σ0, cbf0, kd0, formula="gaussian", σ=0.1);
```

This smooth safety filter is now used to construct a CBF for the double-integrator:
```julia
μ = 5.0
h(x) = h0(x[1:2]) - (0.5/μ)*norm(x[3:4] - k0(x[1:2]))^2
cbf = ControlBarrierFunction(h, s -> s);
```

Finally, we construct a QP-based safety filter for the double-integrator:
```julia
kd(x) = -x[1:2] - 2*x[3:4] # Desired controller for full-order model
k = ReluSafetyFilter(Σ, cbf, kd);
```

If we want, we can then simulate the system under this controller:
```julia
x0 = [-2.1, 2.0, 0.0, 0.0] # Initial condition
T = 15.0 # Length of simulation
sol = simulate(Σ, k, x0, T) # Resulting solution
```

## A Note on Plotting
I make all of my plots using [PGFPlotsX.jl](https://github.com/KristofferC/PGFPlotsX.jl), which is essentially a Julia wrapper for the [PGFPlots LaTeX package](https://www.overleaf.com/learn/latex/Pgfplots_package). The examples in this repo make use of this package, which requires you to have a working LaTeX installation. If you do want to use PGFPlots for plotting, please have a look at the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

## Related Works

### Smooth Safety Filters
- P. Ong and J. Cortes, “[Universal formula for smooth safe stabilization](https://ieeexplore.ieee.org/abstract/document/9030225),” in *Proceeding of the IEEE Conference on Decision and Control*, pp. 2373–2378, 2019.

- M. H. Cohen, P. Ong, G. Bahati, and A. D. Ames, “[Characterizing smooth safety filters via the implicit function theorem](https://ieeexplore.ieee.org/abstract/document/10352951),” *IEEE Control Systems Letters*, 2023.

### Backstepping
- A. J. Taylor, P. Ong, T. G. Molnar, and A. D. Ames, “[Safe backstepping with control barrier functions](https://ieeexplore.ieee.org/abstract/document/9992763),” in *Proceeding of the IEEE Conference on Decision and Control*, pp. 5775–5782, 2022.

### CBFs via Reduced-Order Models
- T. G. Molnar, R. K. Cosner, A. W. Singletary, W. Ubellacker,
and A. D. Ames, “[Model-free safety-critical control for robotic
systems](https://ieeexplore.ieee.org/abstract/document/9652122),” *IEEE Robotics and Automation Letters*, vol. 7, no. 2,
pp. 944–951, 2022.

- T. G. Molnar and A. D. Ames, “[Safety-critical control with
bounded inputs via reduced order models](https://arxiv.org/abs/2303.03247),” in *Proc. Amer. Control Conf.*, pp. 1414–1421, 2023.

- A. W. Singletary, K. Klingebiel, J. Bourne, A. Browning,
P. Tokumaru, and A. D. Ames, “[Comparative analysis of control barrier functions and artificial potential fields for obstacle avoidance](https://ieeexplore.ieee.org/abstract/document/9636670),” in *Proceedings of the IEEE/RSJ International Conference on Intelligent Robots and Systems*, pp. 8129–8136,
2021

- A. Singletary, S. Kolathaya, and A. D. Ames, “[Safety-critical
kinematic control of robotic systems](https://ieeexplore.ieee.org/abstract/document/9319250),” *IEEE Control Systems Letters*,
vol. 6, pp. 139–144, 2022.

- A. W. Singletary, W. Guffey, T. Molnar, R. Sinnet, and
A. D. Ames, “[Safety-critical manipulation for collision-free food
preparation](https://ieeexplore.ieee.org/abstract/document/9834089),” *IEEE Robotics and Automation Letters*, vol. 7,
no. 4, pp. 10954–10961, 2022.