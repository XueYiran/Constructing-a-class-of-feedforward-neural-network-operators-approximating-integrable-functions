import numpy as np
import math


def f_prime(x, A, lambda_, q):
    term = A ** (lambda_ * x) + q ** 2 * A ** (-lambda_ * x) + 2 * q
    return q * lambda_ * math.log(A) / term


def f_double_prime(x, A, lambda_, q):
    term = A ** (lambda_ * x) + q ** 2 * A ** (-lambda_ * x) + 2 * q
    term_squared = term ** 2
    numerator = q * lambda_ ** 2 *math.log(A) ** 2 *(q ** 2 * A ** (-lambda_ * x) - A ** (lambda_ * x))
    return numerator / term_squared


def y_prime(x, A, lambda_, q):
    return 0.5 * (f_prime(x + 1, A, lambda_, q) - f_prime(x - 1, A, lambda_, q))


def y_double_prime(x, A, lambda_, q):
    return 0.5 * (f_double_prime(x + 1, A, lambda_, q) - f_double_prime(x - 1, A, lambda_, q))


def curvature(x, A, lambda_, q):
    y_p = y_prime(x, A, lambda_, q)
    y_pp = y_double_prime(x, A, lambda_, q)
    return np.abs(y_pp) / (1 + y_p ** 2) ** (3 / 2)


q = 1
A = math.e
lambdas = [0.25,0.50,0.75,1.00,1.25,1.50]

for lambda_ in lambdas:
    K = curvature(0, A, lambda_, q)
    print(f"When λ = {lambda_:.2f}, curvature at x = 0 is: {K:.6f}")