import numpy as np
from scipy.integrate import quad
import math

def sigma(x, q,lambda_, A):
    return 1 / (1 + q * A ** (-lambda_ * x))

def phi_sigma(x, q,lambda_, A):
    return 0.5 * (sigma(x + 1, q,lambda_, A) - sigma(x - 1, q,lambda_, A))

def f(x):
    return 1 / (x ** 2 + 1)

def f_0(x):
    if -1 <= x <= 1:
        return f(x)
    elif -2 <= x < -1:
        return -2 * f(-3 - 2 * x) + 3 * f(-2 - x)
    elif 1 < x <= 2:
        return -2 * f(3 - 2 * x) + 3 * f(2 - x)
    else:
        return 0

def K(n, x, q,lambda_, A):
    sum_k = 0
    for k in range(-2 * n, 2 * n):
        x_k = k / n
        x_k_plus_1 = (k + 1) / n
        integral = quad(lambda t: f_0(t), x_k, x_k_plus_1)[0]
        sum_k += integral * phi_sigma(n * x - k, q,lambda_, A)
    return n * sum_k

def print_max_error(ns, As):
    for n in ns:
        for A in As:
            max_error = max(abs(K(n, x, q,lambda_, A) - f(x)) for x in np.linspace(-1, 1, 100))
            print(f"n={n}, A={A:.2f}: Max Error = {max_error:.6f}")

lambda_ = 1
q = 1
ns = [40,50,60,70,80]
As = [1.50, 2.00, 2.50, 3.00, 3.50, 4.00]

print_max_error(ns, As)