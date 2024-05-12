using CairoMakie

kronecker(n::Integer)::Real = n == 0 ? 1 : 0


x = range(-2, step=1, stop=2)
y = [kronecker(t) for t in x]

#Plot
scatter(x,y)