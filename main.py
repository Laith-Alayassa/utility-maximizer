import sympy as sym


def main():
    x, y, l = sym.symbols("x y l")
    utility = (x**0.3) * (y**0.7)
    budget = (1.1 * x) + (y) - 231

    lagrange = setup_lagrange_function(l, utility, budget)
    dx, dy, dl = find_partial_derivs(x, y, l, lagrange)
    deriv_ratio = find_derivative_ratio_equation(dx, dy, l)

    x_in_terms_of_y = get_x_in_terms_of_y(x, deriv_ratio)
    budget_in_terms_of_x = find_budget_equation_in_terms_of_x(x, budget, x_in_terms_of_y)

    y_value = find_y(budget_in_terms_of_x)
    x_value = find_x(y, x_in_terms_of_y, y_value)

    print(f"""
    Utility function = {utility}
    budget constraint = {budget}

    lagrange = {lagrange}

    First order derivatives:

    ∂l/∂x = {dx}

    ∂l/∂y = {dy}

    ∂l/∂λ = {dl}


    X = {x_in_terms_of_y}

    plug x into budget:
    {budget} ==> {budget_in_terms_of_x}

    x = {round(float(x_value), 3)}
    y = {round(float(y_value), 3)}
    """)


def find_x(y, x_in_terms_of_y, y_value):
    return x_in_terms_of_y.replace(y, y_value)


def find_y(budget_in_terms_of_x):
    return sym.solve(budget_in_terms_of_x)[0]


def find_budget_equation_in_terms_of_x(x, budget, x_in_terms_of_y):
    return budget.replace(x, x_in_terms_of_y)


def get_x_in_terms_of_y(x, deriv_ratio):
    return sym.solve(deriv_ratio)[0][x]


def setup_lagrange_function(l, utility, budget):
    lagrange = sym.sympify(utility - l * (budget))
    return lagrange


def find_partial_derivs(x, y, l, lagrange):
    dx, dy, dl = sym.diff(lagrange, x), sym.diff(lagrange, y), sym.diff(lagrange, l)
    return dx, dy, dl


def find_derivative_ratio_equation(dx, dy, l):
    """ because the equation is equal to 0, and I want to take lambda to the right side for both equations
    I will add two lambdas to the right side (one to cancel one from the left side) and the negative is to invert the sign
    since I am moving it to the other side. This will give me something similar to the by-hand solution where lambdas are on the right side.
    """
    x_side_lambda = -2 * dx.coeff(l)
    y_side_lambda = -2 * dy.coeff(l)
    deriv_ratio = sym.Eq(dx / dy, x_side_lambda/y_side_lambda)
    return deriv_ratio


if __name__ == '__main__':
    main()
