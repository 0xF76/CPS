using CairoMakie

x = range(0.25, step=0.001, length=256)
y = [2 * sin(2 * π * 25 * t) for t in x]

# Plot
scatter(x, y, color=:blue)