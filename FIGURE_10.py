import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math

def binomial_coefficient(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))

def m_r(x, r):
    if x <= -r / 2:
        return 0
    elif x >= r / 2:
        return 1
    else:
        segment = int(np.floor(x + r / 2))
        result = (x + r / 2) ** r
        for k in range(1, segment + 1):
            if k % 2 == 0:
                result += binomial_coefficient(r, k) * (x + r / 2 - k) ** r
            else:
                result -= binomial_coefficient(r, k) * (x + r / 2 - k) ** r
        return result / math.factorial(r)

def phi_m_r(x, r):
    return 0.5 * (m_r(x + 1, r) - m_r(x - 1, r))

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

def K(n, x, r):
    sum_k = 0
    for k in range(-2 * n, 2 * n):
        x_k = k / n
        x_k_plus_1 = (k + 1) / n
        integral = quad(lambda t: f_0(t), x_k, x_k_plus_1)[0]
        sum_k += integral * phi_m_r(n * x - k, r)
    return n * sum_k

def plot_max_error(ns,rs):
    max_errors = {}
    for r in rs:
        errors = []
        for n in ns:
            def integrand(x):
                return abs(K(n, x, r) - f(x))
            integral, _ = quad(integrand, -1, 1, limit=100, epsabs=1e-6, epsrel=1e-6)
            error = integral
            errors.append(error)
        max_errors[r] = errors

    plt.figure(figsize=(10, 6))
    for r, errors in max_errors.items():
        plt.plot(ns, errors, label=f'r={r}', marker='o', linestyle='-')
    plt.xlabel('n')
    plt.ylabel('Max Error')
    plt.legend()
    plt.grid(False)
    plt.xticks(ns)
    plt.show()

ns = [40,50,60,70,80]
rs = [2, 3, 4, 5, 6, 7]
plot_max_error(ns, rs)