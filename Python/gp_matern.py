import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel

# Load data file
data = np.loadtxt("output/shannon_geneorientation_.dat", comments="!")

H    = data[:,0]   # Shannon entropy
frac = data[:,1]   # fraction convergent/divergent
std  = data[:,2]   # spacing std

# input coordinates for GP
X = np.column_stack((std, frac))
y = H

# Define Matérn kernel
kernel = (
    ConstantKernel(1.0, (1e-3, 1e3))
    * Matern(length_scale=[20, 0.3],
             length_scale_bounds=(1e-2, 1e3),
             nu=1.5)
    + WhiteKernel(noise_level=1e-5,
                  noise_level_bounds=(1e-8, 1e-1))
)

# Fit Gaussian Process
gp = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=10,
        normalize_y=True
)
gp.fit(X, y)

print("Optimized kernel:", gp.kernel_)

# Create grid
xi = np.linspace(std.min(), std.max(), 200)
yi = np.linspace(frac.min(), frac.max(), 200)
XI, YI = np.meshgrid(xi, yi)

Xgrid = np.column_stack([XI.ravel(), YI.ravel()])

# Gaussian Process prediction
ZI, sigma = gp.predict(Xgrid, return_std=True)
ZI = ZI.reshape(XI.shape)
sigma = sigma.reshape(XI.shape)


# Find location of Shannon entropy minimum
min_index = np.argmin(ZI)
min_std = Xgrid[min_index, 0]
min_frac = Xgrid[min_index, 1]
min_H = ZI.ravel()[min_index]

print("\nPredicted minimum entropy:")
print("std =", min_std)
print("fraction =", min_frac)
print("H =", min_H)

# Plot heatmap
plt.figure()

contour = plt.contourf(XI, YI, ZI, 50)

plt.scatter(std, frac, c=H, edgecolor='k',)

plt.scatter(min_std, min_frac,
            marker="*", s=200)

plt.colorbar(contour, label="Shannon entropy H")

plt.xlabel("Std of inter-gene distance")
plt.ylabel("Fraction of convergent/divergent neighbours")
plt.title("Gaussian Process entropy landscape (Matern kernel)")

plt.show()