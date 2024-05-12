function running_mean(x::AbstractVector, M::Integer)::Vector
    rm = []
    for i in 1:length(x)
        sum = 0
        for j in (i-M):(i+M)
            if j < 1 || j > length(x)
                sum += 0
            else
                sum += x[j]
            end 
        end
        push!(rm, sum/(2*M+1))
    end 
    return rm
end


signal = 1:10
running_mean(signal, 2)