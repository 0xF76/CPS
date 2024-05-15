using CairoMakie

cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T/2 ? 1 : 0


x = range(-2, step=0.01, stop=2)
y = [cw_rectangular(t) for t in x]

#Plot
lines(x, y, color=:blue)