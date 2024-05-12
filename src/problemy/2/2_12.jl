using CairoMakie

square_wave(t::Real)::Real =  sign(sin(2*π*t))

x = range(-2, step=0.01, stop=2)
y = [square_wave(t) for t in x]

#Plot
lines(x, y, color=:blue)