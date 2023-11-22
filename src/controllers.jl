abstract type Controller end
abstract type FeedbackController <: Controller end
abstract type SafetyFilter <: FeedbackController end
abstract type QPController <: FeedbackController end
abstract type MPController <: FeedbackController end