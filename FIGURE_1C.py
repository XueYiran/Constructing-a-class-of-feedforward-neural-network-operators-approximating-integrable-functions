import numpy as np
import matplotlib.pyplot as plt


def sigma(x, q):
    return 1 / (1 + q * np.exp(-x))

def phi_sigma(x, q):
    return 0.5 * (sigma(x + 1, q) - sigma(x - 1, q))

x = np.linspace(-7, 7, 700)

plt.figure(figsize=(10, 6))
ax = plt.gca()

qs = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75]
colors = ['r', 'g', 'b', 'm', 'c', 'y']

for i, q in enumerate(qs):
    y = phi_sigma(x, q)
    plt.plot(x, y, label=f'q = {q:.2f}', color=colors[i])

x_labels = np.arange(-7, 8, 1)
plt.xticks(x_labels)


plt.yticks(np.arange(0.0, 0.26, 0.05))
plt.ylim(0.0, 0.25)

plt.legend()

plt.grid(False)

plt.show()