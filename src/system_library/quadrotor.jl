"""
    QuadrotorControlAffine <: ControlAffineSystem

Control affine representation of a 3D quadrotor.
"""
struct QuadrotorControlAffine <: ControlAffineSystem
    n::Int
    m::Int
    f::Function
    g::Function
end
Quadrotor() = QuadrotorControlAffine()
Quadrotor(::Type{ControlAffineSystem}) = QuadrotorControlAffine()

"""
    QuadrotorControlAffine()

Construct Quadrotor as a control affine system. Dynamics based on MATLAB code provided by Ryan Cosner.
"""
function QuadrotorControlAffine()
    # Constants
    mass = 1.0 # Drone mass
    grav = 9.81 # Acceleration due to gravity
    ez = [0.0, 0.0, 1.0] # Unit vector in z direction
    unit_z_body_quat = quat(0.0, 0.0, 0.0, 1.0)

    # System dimensions
    """
        state = [   
            1, x , world-frame position
            2, y , world-frame position
            3, z , world-frame position
            4, qw  , world-frame body quaternion
            5, qx  , world-frame body quaternion
            6, qy  , world-frame body quaternion
            7, qz  , world-frame body quaternion 
            8, vx  , world-frame body velocity
            9, vy  , world-frame body velocity
            10, vz  , world-frame body velocity 
        ]
    """
    n = 10 
    m = 4

    # Drift dynamics
    f(x) = [x[8:10]; zeros(4); -grav*ez]

    # Control directions
    function g(x)
        ctrl_dir = zeros(n, m)
        quat_world = normalize(quat(x[4:7]...))
        thrust_dir = normalize(vectorize3tuple(imag_part(quat_world*unit_z_body_quat*inv(quat_world))))
        world_qw, world_qx, world_qy, world_qz = parts(quat_world)
        angle_rate_input_matrix = [
                                    0.0  world_qx  world_qy  world_qz;
                                    0.0  world_qw -world_qz  world_qy;
                                    0.0  world_qz  world_qw -world_qx;
                                    0.0  world_qy  world_qx  world_qw
                                ]
        
        # Allocate control directions
        ctrl_dir[8:10, 1] = thrust_dir
        ctrl_dir[4:7, :] = angle_rate_input_matrix/mass

        return ctrl_dir
    end

    return QuadrotorControlAffine(n, m, f, g)
end

"""
    dynamics(Σ::QuadrotorControlAffine, x)
    dynamics(Σ::QuadrotorControlAffine, x, u)

Compute `ẋ` for the control affine representation of a 3D quadrotor.
"""
function dynamics(Σ::QuadrotorControlAffine, x)
    return Σ.f(x)
end
function dynamics(Σ::QuadrotorControlAffine, x, u)
    q = quat(x[4:7]...)
    R = rotmatrix_from_quat(quat(q))
    u[2:4] = R'*u[2:4]
    ẋ = Σ.f(x) + Σ.g(x)*u
    
    return ẋ
end

"""
    simulate(Σ::QuadrotorControlAffine, x, T)
    simulate(Σ::QuadrotorControlAffine, x, k::FeedbackController, T)

Simulate quadrotor with callback to terminate sim if we hit the ground.
"""
function simulate(Σ::QuadrotorControlAffine, x, T)
    # Create callback to terminate sim when we hit ground
    condition(u, t, integrator) = u[3]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    prob = ODEProblem((x,p,t) -> Σ.f(x), x, (0,T))
    sol = solve(prob, callback=cb)

    return sol
end

function simulate(Σ::QuadrotorControlAffine, x, k::FeedbackController, T)
    # Construct dynamics function
    f(x) = dynamics(Σ, x, k(x))

    # Create callback to terminate sim when we hit ground
    condition(u, t, integrator) = u[3]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    prob = ODEProblem((x,p,t) -> f(x), x, (0,T))
    sol = solve(prob, callback=cb)

    return sol
end

"""
    rotmatrix_from_quat(q::Quaternion)

Convert a quaternion into a rotation matrix. Taken directly from: 
https://juliageometry.github.io/Quaternions.jl/stable/examples/rotations/
"""
function rotmatrix_from_quat(q::Quaternion)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    r = [1 - (yy + zz)     xy - sz     xz + sy;
            xy + sz   1 - (xx + zz)    yz - sx;
            xz - sy      yz + sx  1 - (xx + yy)]
    return r
end

function rotmatrix_from_quat(s::Float64, v1::Float64, v2::Float64, v3::Float64)
    q = quat(s, v1, v2, v3)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    r = [1 - (yy + zz)     xy - sz     xz + sy;
            xy + sz   1 - (xx + zz)    yz - sx;
            xz - sy      yz + sx  1 - (xx + yy)]
    return r
end

"""
    parts(q::Quaternion)

Extract quaternion parts.
"""
parts(q::Quaternion) = [real(q); imag_part(q)...]

"""
    vectorize3tuple(v)

Convert a tuple of 3 elements into a vector.
"""
vectorize3tuple(v) = [v[1], v[2], v[3]]

"""
    vectorize_skew_matrix(M)

Convert skew-symmetric matrix into a 3 vector.
"""
vectorize_skew_matrix(M) = [M[3,2], M[1,3], M[3,1]]

"""
    DiffFlatQuadController <: FeedbackController

Nominal controller used for quadrotor based on differential flatness.
"""
struct DiffFlatQuadController <: FeedbackController
    get_input::Function
end
(k::DiffFlatQuadController)(x) = k.get_input(x)

"""
    DiffFlatQuadController()

Construct for quad's differential flatness controller.
"""
function DiffFlatQuadController(; xd=[2.0, 1.0, 3.0])
    # Useful stuff to have
    mass = 1.0 # Drone mass
    grav = 9.81 # Acceleration due to gravity
    ez = [0.0, 0.0, 1.0] # Unit vector in z direction

    # Desired position and velocity; TODO: make these parameters of controller?
    # xd = [2.0, 1.0, 3.0] # Desired position
    vd = zeros(3) # Desired velocity

    # Controller gains; TODO: make these parameters of controller?
    Kx = 1 
    Kv = 1 
    Kangle = 1 

    # Control function
    function get_input(x)
        # Get orientation of quadrotor
        q = quat(x[4:7]...)
        R = rotmatrix_from_quat(quat(q))

        # Position and velocity error
        ex = x[1:3] - xd
        ev = x[8:10] - vd

        b1_des = [1.0, 0.0, 0.0]
        b3_des = normalize(-Kx*ex - Kv*ev  + mass*grav*ez)
        b2_des = normalize(cross(b3_des, b1_des))

        R_des = hcat(cross(b2_des, b3_des), b2_des, b3_des)
        angle_error = 1.0 / 2.0 * vectorize_skew_matrix(R_des'*R - R'*R_des)
        angle_rate_des = - Kangle * angle_error

        thrust_force = dot(-Kx*ex - Kv*ev + mass*grav*ez, R*ez)

        u = [thrust_force; angle_rate_des]

        return u
    end

    return DiffFlatQuadController(get_input)
end