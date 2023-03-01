import sympy as sym


# px * x + py * y = i
# a * x  + b  * y - i = 0


# variables
a, x, b, y, i, l = sym.symbols("a x b y i l")
utility = (x ** 0.5) * (y ** 0.5)
budget = (1 * x) + (0.5 * y - 2)

# set up lagrange
lagrange = sym.sympify(utility - l * (budget))

# deriv
dx = sym.diff(lagrange, x)
dy = sym.diff(lagrange, y)
dl = sym.diff(lagrange, l)


# ratios
x_side_lambda = -2 * dx.coeff(l)  # takes lambda to the other side
y_side_lambda = -2 * dy.coeff(l)

deriv_ratio = sym.Eq(dx / dy, x_side_lambda/y_side_lambda)

X_in_terms_of_Y = sym.solve(deriv_ratio)[0][x]

budget_in_terms_of_x = budget.replace(x, X_in_terms_of_Y)

# print(budget_in_terms_of_x)

y_value = sym.solve(budget_in_terms_of_x)[0]

x_value = X_in_terms_of_Y.replace(y, y_value)


print(f"""
      Utility function = {utility}
      budget constraint = {budget}
      
      lagrange = {lagrange}
      
      First order derivatives:
      
      ∂l/∂x = {dx}
      
      ∂l/∂y = {dx}
      
      ∂l/∂λ = {dl}
      
      
      ratio = {dx/dy} == {x_side_lambda}λ/{y_side_lambda}λ
      
      X = {X_in_terms_of_Y}
      
      
      plug x into budget:
      {budget} ==> {budget_in_terms_of_x}
      
      y = {y_value}
      x = {x_value}
      """)
