abstract type Controller end
abstract type FeedbackController <: Controller end
abstract type SafetyFilter <: FeedbackController end
abstract type QPController <: FeedbackController end

"""
    StateFeedbackController <: FeedbackController

Custom state feedback controller. 
"""
struct StateFeedbackController <: FeedbackController
    get_input::Function
end
(k::StateFeedbackController)(x) = k.get_input(x)
(k::StateFeedbackController)(q, q̇) = k.get_input(q, q̇)