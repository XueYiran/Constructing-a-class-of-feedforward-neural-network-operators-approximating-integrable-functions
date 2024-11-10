import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math

def sigma(x, q, lambda_, A):
    return 1 / (1 + q * A ** (-lambda_ * x))

def phi_sigma(x, q, lambda_, A):
    return 0.5 * (sigma(x + 1, q, lambda_, A) - sigma(x - 1, q, lambda_, A))

def f(x):
    if -1 <= x < 0:
        return 1 / (x ** 2 + 1)
    elif 0 <= x <= 1:
        return math.sin(x)

def f_0(x):
    if -1 <= x <= 1:
        return f(x)
    elif -2 <= x < -1:
        return -2 * f(-3 - 2 * x) + 3 * f(-2 - x)
    elif 1 < x <= 2:
        return -2 * f(3 - 2 * x) + 3 * f(2 - x)
    else:
        return 0


def K(n, x, q, lambda_, A):
    sum_k = 0
    for k in range(-2 * n, 2 * n):
        x_k = k / n
        x_k_plus_1 = (k + 1) / n
        integral = quad(lambda t: f_0(t), x_k, x_k_plus_1)[0]
        sum_k += integral * phi_sigma(n * x - k, q, lambda_, A)
    return n * sum_k


def plot_max_error(q, A, ns, lambda_s):
    max_errors = {}
    for lambda_ in lambda_s:
        errors = []
        for n in ns:
            def integrand(x):
                return abs(K(n, x, q, lambda_, A) - f(x))
            integral, _ = quad(integrand, -1, 1, limit=100, epsabs=1e-6, epsrel=1e-6)
            error = integral
            errors.append(error)
        max_errors[lambda_] = errors

    plt.figure(figsize=(10, 6))
    for lambda_, errors in max_errors.items():
        plt.plot(ns, errors, label=f'λ={lambda_:.2f}', marker='o', linestyle='-')
    plt.xlabel('n')
    plt.ylabel('Max Error')
    plt.legend()
    plt.grid(False)
    plt.xticks(ns)
    plt.show()

q = 1
A = math.e
ns = [40,50,60,70,80]
lambda_s = [0.25,0.50,0.75,1.00,1.25,1.50]
plot_max_error(q, A, ns, lambda_s)