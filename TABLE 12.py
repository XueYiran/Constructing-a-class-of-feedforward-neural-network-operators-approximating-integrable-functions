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

ns = [40,50,60,70,80]
rs = [2, 3, 4, 5, 6, 7]

for n in ns:
    for r in rs:
        max_error = max(abs(K(n, x, r) - f(x)) for x in np.linspace(-1, 1, 100))
        print(f"n={n}, r={r}: Max Error = {max_error:.7f}")




