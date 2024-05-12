using CairoMakie

pulse_wave(t::Real, ρ::Real=0.2)::Real = (t - floor(t)) ≤ ρ ? 1 : 0

x = range(-2, step=0.01, stop=2)
y = [pulse_wave(t) for t in x]

#Plot
lines(x, y, color=:blue)
