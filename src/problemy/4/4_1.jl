
using CairoMakie

t = range(0,step = 1/250, stop = 0.1)
x = [sin(200 * π * i) for i in t]


scatter(t,x)


