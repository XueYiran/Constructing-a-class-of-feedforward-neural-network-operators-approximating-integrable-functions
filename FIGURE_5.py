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

def phi_m_r(x, r):
    return 0.5 * (m_r(x + 1, r) - m_r(x - 1, r))

x = np.linspace(-4, 4, 400)
plt.figure(figsize=(10, 6))

for r in range(2, 8):
    y = [phi_m_r(xi, r) for xi in x]
    plt.plot(x, y, label=f'r={r}')


plt.legend()
plt.grid(False)
plt.show()