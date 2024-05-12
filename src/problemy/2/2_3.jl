white_noise(n) = √0.25 * randn(n)


x = white_noise(100)

for i in 1:1000
    print(x[i], ",")
end




