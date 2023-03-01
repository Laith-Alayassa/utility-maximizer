import unittest
from main import solve_equation
import sympy as sym


class TestSum(unittest.TestCase):
    def test_solve_equaiton(self):
        self.assertEqual(solve_equation(sym.sympify("x+3"), sym.Symbol("x")), [-3])


if __name__ == '__main__':
    unittest.main()
