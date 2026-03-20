import numpy as np
import matplotlib.pyplot as plt
import os

# Create results folder
os.makedirs("results", exist_ok=True)

# Grid
N = 100
dx = 1.0
dt = 0.01
steps = 200

# Parameters
Dc, Dn, Ds = 0.1, 0.01, 0.1
gamma = 0.5
alpha = 0.05
beta = 0.01
mu0 = 0.3
K = 0.1

# Initialize fields
c = np.ones((N, N)) * 1.0
n = np.zeros((N, N))
s = np.zeros((N, N))

# Initial condition (small bacteria seed)
n[N//2-1:N//2+1, N//2-1:N//2+1] = 1.0


def laplacian(Z):
    return (
        np.roll(Z, 1, axis=0) + np.roll(Z, -1, axis=0) +
        np.roll(Z, 1, axis=1) + np.roll(Z, -1, axis=1) -
        4 * Z
    ) / dx**2


# Store growth
colony_size = []

# =========================
# MAIN SIMULATION
# =========================
for t in range(steps):
    mu = mu0 * c / (K + c)

    c_new = c + dt * (Dc * laplacian(c) - gamma * n * c)
    n_new = n + dt * (Dn * laplacian(n) + mu * n)
    s_new = s + dt * (Ds * laplacian(s) + alpha * n - beta * s)

    # Prevent negative / explosion
    c = np.clip(c_new, 0, 2)
    n = np.clip(n_new, 0, 10)
    s = np.clip(s_new, 0, 5)

    colony_size.append(np.sum(n))


# =========================
# 1. COLONY DENSITY PLOT
# =========================
plt.figure()
plt.imshow(n, origin="lower")
plt.title("Bacterial Density")
plt.colorbar()
plt.tight_layout()
plt.savefig("results/colony_growth.png")
plt.show()


# =========================
# 2. PARAMETER HEATMAP
# =========================
Dn_vals = np.linspace(0.005, 0.02, 6)
Dc_vals = np.linspace(0.05, 0.25, 6)

heatmap = np.zeros((len(Dc_vals), len(Dn_vals)))

for i, Dc_val in enumerate(Dc_vals):
    for j, Dn_val in enumerate(Dn_vals):

        c_temp = np.ones((N, N))
        n_temp = np.zeros((N, N))
        n_temp[N//2, N//2] = 1.0

        for _ in range(100):
            mu = mu0 * c_temp / (K + c_temp)

            c_temp = np.clip(c_temp + dt * (Dc_val * laplacian(c_temp) - gamma * n_temp * c_temp), 0, 2)
            n_temp = np.clip(n_temp + dt * (Dn_val * laplacian(n_temp) + mu * n_temp), 0, 10)

        heatmap[i, j] = np.sum(n_temp)

plt.figure()
plt.imshow(heatmap, origin="lower", aspect='auto')
plt.colorbar(label="Total Bacteria")
plt.xticks(range(len(Dn_vals)), np.round(Dn_vals, 3))
plt.yticks(range(len(Dc_vals)), np.round(Dc_vals, 3))
plt.xlabel("Dn (Motility)")
plt.ylabel("Dc (Diffusion)")
plt.title("Parameter Sensitivity Heatmap")
plt.tight_layout()
plt.savefig("results/heatmap.png")
plt.show()


# =========================
# 3. GROWTH CURVE
# =========================
plt.figure()
plt.plot(colony_size)
plt.xlabel("Time Step")
plt.ylabel("Total Bacteria")
plt.title("Colony Growth Over Time")
plt.grid()
plt.tight_layout()
plt.savefig("results/growth_curve.png")
plt.show()
