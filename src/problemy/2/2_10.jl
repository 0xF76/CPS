using CairoMakie

sawtooth_wave(t::Real)::Real = -(t-floor(t))


x = range(-2, step=0.01, stop=2)
y = [sawtooth_wave(t) for t in x]

#Plot
lines(x, y, color=:blue)