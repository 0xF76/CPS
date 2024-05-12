peak2peak(x::AbstractVector)::Real = maximum(x) - minimum(x)


signal = 1:10
peak2peak(signal)