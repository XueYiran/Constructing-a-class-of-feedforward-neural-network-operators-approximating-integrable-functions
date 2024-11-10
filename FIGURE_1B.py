import numpy as np
import matplotlib.pyplot as plt


def sigma(x, A):
    return 1 / (1 + A ** (-x))

def phi_sigma(x, A):
    return 0.5 * (sigma(x + 1, A) - sigma(x - 1, A))

x = np.linspace(-15, 15, 1500)

plt.figure(figsize=(10, 6))
ax = plt.gca()

As = [1.50, 2.00, 2.50, 3.00, 3.50, 4.00]
colors = ['r', 'g', 'b', 'm', 'c', 'y']

for i, A in enumerate(As):
    y = phi_sigma(x, A)
    plt.plot(x, y, label=f'A = {A:.2f}', color=colors[i])

x_labels = np.arange(-15, 16, 3)
plt.xticks(x_labels)

plt.yticks(np.arange(0.0, 0.36, 0.05))
plt.ylim(0.0, 0.35)

plt.legend()

plt.grid(False)

plt.show()