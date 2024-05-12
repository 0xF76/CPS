using CairoMakie

x = range(5, step=1/2048, stop = 10)
y = [(0.25 * exp(im*π/2*t + π)) for t in x]

real_y = real.(y)
imag_y = imag.(y)

scatter(x, real_y, color=:blue, markersize=3, label="real")
scatter!(x, imag_y, color=:red, markersize=3, label="imaginary")
axislegend()

current_figure()