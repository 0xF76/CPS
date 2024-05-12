energy(x::AbstractVector)::Real = sum(x.^2)

signal = 1:10

energy(signal)