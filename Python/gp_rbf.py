import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from scipy.optimize import minimize

# Load file with shannon entropy for specific spacing and fraction of convergent/divergent genes
data = np.loadtxt("output/shannon_geneorientation_.dat", comments="!")

# File structure
H    = data[:,0]
frac = data[:,1]
std  = data[:,2]

# 2D vector of spacing and fraction of convergent/divergent genes
X = np.column_stack([std, frac])


# Define GP model
kernel = (
    ConstantKernel(1.0, (1e-3, 1e3))
    * RBF(length_scale=[20, 0.3], length_scale_bounds=(1e-2, 1e2))
    + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-8, 1e-1))
)

gp = GaussianProcessRegressor(
    kernel=kernel,
    n_restarts_optimizer=10,
    normalize_y=True
)

gp.fit(X, H)

print("Learned kernel:")
Print(gp.kernel_)

# Find minimum of GP parameter space
def gp_function(x):
    return gp.predict(np.array(x).reshape(1,-1))[0]

x0 = [std.mean(), frac.mean()]

res = minimize(
    gp_function,
    x0,
    bounds=[
        (std.min(), std.max()),
        (frac.min(), frac.max())
    ]
)

print("\nEstimated minimum:")
print("STD =", res.x[0])
print("FRAC =", res.x[1])
print("H =", res.fun)


# Create smooth grid for plot
xi = np.linspace(std.min(), std.max(), 200)
yi = np.linspace(frac.min(), frac.max(), 200)

XI, YI = np.meshgrid(xi, yi)
grid_points = np.column_stack([XI.ravel(), YI.ravel()])

ZI = gp.predict(grid_points).reshape(XI.shape)

# Plot smooth Shannon entropy landscape
plt.figure(figsize=(7,5))

contour = plt.contourf(XI, YI, ZI, 100)
plt.colorbar(contour, label="Shannon entropy H")

# data points
plt.scatter(std, frac, c=H, edgecolor='k', s=60)

# minimum point
plt.scatter(res.x[0], res.x[1], marker='*', s=200)

plt.xlabel("Std of inter-gene distance")
plt.ylabel("Fraction of convergent/divergent neighbours")
plt.title("Gaussian Process entropy landscape - RBF")

plt.show()