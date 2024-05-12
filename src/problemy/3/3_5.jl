energy(x::AbstractVector)::Real = sum(x.^2)
power(x::AbstractVector)::Real = energy(x)/length(x)
rms(x::AbstractVector)::Real = √power(x)

signal = 1:10

rms(signal)