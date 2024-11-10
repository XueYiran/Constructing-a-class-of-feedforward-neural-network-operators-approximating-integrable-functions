import numpy as np
import matplotlib.pyplot as plt
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
            if k%2==0:
                result += binomial_coefficient(r, k) * (x + r / 2 - k) ** r
            else:
                result -= binomial_coefficient(r, k) * (x + r / 2 - k) ** r
        return result / math.factorial(r)

def derivative1_m_r(x, r, h=1e-5):
    return (m_r(x + h, r) - m_r(x - h, r)) / (2 * h)
def derivative2_m_r(x,r, h=1e-5):
    df_x_plus_h = derivative1_m_r(x + h / 2, r)
    df_x_minus_h = derivative1_m_r(x - h / 2, r)
    return (df_x_plus_h - df_x_minus_h) / h

def f_prime(x, r):
    return derivative1_m_r(x, r, h=1e-5)

def f_double_prime(x, r):
    return derivative2_m_r(x,r, h=1e-5)

def y_prime(x, r):
    return 0.5 * (f_prime(x + 1, r) - f_prime(x - 1, r))

def y_double_prime(x, r):
    return 0.5 * (f_double_prime(x + 1, r) - f_double_prime(x - 1, r))

def curvature(x, r):
    y_p = y_prime(x,r)
    y_pp = y_double_prime(x, r)
    return np.abs(y_pp) / (1 + y_p ** 2) ** (3 / 2)

rs=[2,3,4,5,6,7]

for r in rs:
    K = curvature(0, r)
    print(f"When r = {r}, curvature at x = 0 is: {K:.6f}")
