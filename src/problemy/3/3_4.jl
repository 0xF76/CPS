energy(x::AbstractVector)::Real = sum(x.^2)
power(x::AbstractVector)::Real = energy(x)/length(x)

signal = 1:10

power(signal)