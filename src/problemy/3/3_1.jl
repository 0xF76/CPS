mean(x::AbstractVector)::Number = sum(x)/length(x)


signal = 1:10
mean(signal)