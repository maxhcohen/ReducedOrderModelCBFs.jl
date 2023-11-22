struct IDSafetyFilter <: FeedbackController
    get_input::Function
end
(k::IDSafetyFilter)(q, q̇) = k.get_input(q, q̇)

function IDSafetyFilter(Σ::RoboticSystem, k::Function, K)
    dk(q) = ForwardDiff.jacobian(k ,q)
    get_input(q, q̇) = inverse_dynamics(Σ, q, q̇, dk(q)*q̇) - K*Σ.M(q)*(q̇ - k(q))

    return IDSafetyFilter(get_input)
end
