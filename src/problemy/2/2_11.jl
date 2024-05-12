using CairoMakie

triangular_wave(t::Real)::Real = 4*abs(t-floor(t+1/2)) - 1


x = range(-2, step=0.01, stop=2)
y = [triangular_wave(t) for t in x]

#Plot
lines(x, y, color=:blue)