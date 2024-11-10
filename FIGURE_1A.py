import numpy as np
import matplotlib.pyplot as plt

def sigma(x, lambda_):
    return 1 / (1 + np.exp(-lambda_ * x))

def phi_sigma(x, lambda_):
    return 0.5 * (sigma(x + 1, lambda_) - sigma(x - 1, lambda_))

x = np.linspace(-15, 16, 1500)

plt.figure(figsize=(10, 6))
ax = plt.gca()

lambdas = [0.25,0.50,0.75,1.00,1.25,1.50]
colors = ['r', 'g', 'b', 'm', 'c', 'y']

for i, lambda_ in enumerate(lambdas):
    y = phi_sigma(x, lambda_)
    plt.plot(x, y, label=f'λ = {lambda_:.2f}', color=colors[i])

x_labels = np.arange(-15, 16, 3)
plt.xticks(x_labels)

plt.yticks(np.arange(0.0, 0.36, 0.05))
plt.ylim(0.0, 0.35)

plt.legend()

plt.grid(False)

plt.show()
