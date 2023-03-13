import sympy as sym


def solve_utility_max():
    global utility_func, budget_line, lagrange_func
    x, y = sym.symbols("x y")

    utility_func = (x**0.5) * (y**0.5)
    budget_line = (1 * x) + (0.5 * y) - 2

    lagrange_func = setup_lagrange_func()
    partial_derivs = find_partial_derivs()
    deriv_ratio = find_deriv_ratio_eq(partial_derivs)

    x_in_terms_of_y = get_x_in_terms_of_y(deriv_ratio)
    budget_in_y = find_budget_eq_in_terms_of_y(x_in_terms_of_y)

    y_value = find_y(budget_in_y)
    x_value = find_x(x_in_terms_of_y, y_value)

    printable_steps = {
        "partial_derivs": partial_derivs,
        "x_in_terms_of_y": x_in_terms_of_y,
        "budget_in_y": budget_in_y,
        "y_value": y_value,
        "x_value": x_value
    }

    solution_text = get_solution_as_string(printable_steps)
    print(get_solution_as_string(printable_steps))
    return solution_text


def find_budget_eq_in_terms_of_y(x_in_terms_of_y):
    x = sym.Symbol("x")
    return budget_line.replace(x, x_in_terms_of_y)


def setup_lagrange_func():
    l = sym.Symbol("l")
    lagrange = sym.sympify(utility_func - l * budget_line)
    return lagrange


def find_partial_derivs():
    x, y, l = sym.symbols("x y l")
    dx, dy, dl = sym.diff(lagrange_func, x), sym.diff(lagrange_func, y), sym.diff(lagrange_func, l)
    return {"dx": dx, "dy": dy, "dl": dl}


def get_x_in_terms_of_y(deriv_ratio):
    x = sym.Symbol("x")
    return sym.solve(deriv_ratio)[0][x]


def find_deriv_ratio_eq(partial_derivs):
    l = sym.Symbol("l")
    dx, dy = partial_derivs["dx"], partial_derivs["dy"]

    """
    because the equation is equal to 0, and I want to take lambda to the right side for both equations
    I will add two lambdas to the right side (one to cancel one from the left side) and the negative is to invert the sign
    since I am moving it to the other side. This will give me something similar to the by-hand solution where lambdas are on the right side.
    """
    x_lambda_coeff = -2 * dx.coeff(l)
    y_lambda_coeff = -2 * dy.coeff(l)

    mrs = dx / dy
    price_ratio = x_lambda_coeff/y_lambda_coeff

    deriv_ratio = sym.Eq(mrs, price_ratio)
    return deriv_ratio


def find_y(budget_in_terms_of_x):
    return sym.solve(budget_in_terms_of_x)[0]


def find_x(x_in_terms_of_y, y_value):
    y = sym.Symbol("y")
    return x_in_terms_of_y.replace(y, y_value)


def get_solution_as_string(steps):
    partial_derivs = steps["partial_derivs"]
    return (
        f"""
    Utility function = {utility_func}  
    budget constraint = {budget_line} 
 
    lagrange = {lagrange_func} 
 
    First order derivatives: 
    ∂l/∂x = {partial_derivs['dx']} 
    ∂l/∂y = {partial_derivs['dy']} 
    ∂l/∂λ = {partial_derivs['dl']} 
 
    Budget Line: {budget_line} = 0   # plug x into budget 
 
    X = {steps["x_in_terms_of_y"]} 
    Budget Line: {steps["budget_in_y"]} = 0 
 
    Solve for Y: 
    Y = {round(float(steps["y_value"]), 3)} 
 
    Plug Y in X equation: 
    X = {steps["x_in_terms_of_y"]} 
    Y = {round(float(steps["y_value"]), 3)} 
    
    ∴ X = {round(float(steps["x_value"]), 3)}   □ 
    """).replace("**", "^")


if __name__ == '__main__':
    solve_utility_max()
