from sympy import symbols, Eq, solve

vert1x, vert1y, vert2x, vert2y, px, py, x, t = symbols("vert1.x vert1.y vert2.x vert2.y projected.x projected.y x t")

eq1 = Eq(px, vert1x + (vert2x - vert1x) * x)
eq2 = Eq(py, vert1y + (vert2y - vert1y) * t)

solution = solve((eq1, eq2), (x, t))

print(solution)
