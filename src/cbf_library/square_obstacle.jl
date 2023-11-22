struct SquareObstacle
    x
    r
end

function ControlBarrierFunction(O::SquareObstacle, α::Function, σ)
    xc = O.x
    r = O.r
    x1 = xc + [0.0, r] # Top of square
    x2 = xc + [r, 0.0] # Right of square
    x3 = xc - [0.0, r] # Bottom of square
    x4 = xc - [r, 0.0] # Left of square
    n1 = [0.0, 1.0] # Normal out of top
    n2 = [1.0, 0.0] # Normal out of right
    n3 = [0.0, -1.0] # Normal out of bottom
    n4 = [-1.0, 0.0] # Normal out of left
    h1(x) = x[2] - x1[2]
    h2(x) = x[1] - x2[1]
    h3(x) = x3[2] - x[2]
    h4(x) = x4[1] - x[1]
    cbf1 = ControlBarrierFunction(h1, α)
    cbf2 = ControlBarrierFunction(h2, α)
    cbf3 = ControlBarrierFunction(h3, α)
    cbf4 = ControlBarrierFunction(h4, α)

    return smooth_disjunction([cbf1, cbf2, cbf3, cbf4], α, σ)
end