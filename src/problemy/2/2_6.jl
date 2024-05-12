using CairoMakie

cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - 1/T*abs(t) : 0

x = range(-2, step=0.01, stop=2)
y = [cw_triangle(t) for t in x]

#Plot
lines(x, y, color=:blue)