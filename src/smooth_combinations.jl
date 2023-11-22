function smooth_conjunction(H::Vector{ControlBarrierFunction}, α, σ)
    h(x) = -σ*log(sum([exp(-cbf.h(x)/σ) for cbf in H]))
    cbf = ControlBarrierFunction(h, α)

    return cbf
end

function smooth_disjunction(H::Vector{ControlBarrierFunction}, α, σ)
    h(x) = σ*log(sum([exp(cbf.h(x)/σ) for cbf in H])) - σ*log(length(H))
    cbf = ControlBarrierFunction(h, α)

    return cbf
end