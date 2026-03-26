import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Load data file
data = np.loadtxt("output/shannon_geneorientation_.dat", comments="!")

H    = data[:,0] # Shannon entropy
std  = data[:,2] # standard deviation in inter-gene spacing
frac = data[:,1] # fraction of convergent/divergent genes

# Create interpolation grid
xi = np.linspace(std.min(), std.max(), 100)
yi = np.linspace(frac.min(), frac.max(), 100)
XI, YI = np.meshgrid(xi, yi)

# Interpolate Shannon entropy parameter space
try:
    ZI = griddata((std, frac), H, (XI, YI), method='linear')
except Exception:
    ZI = griddata((std, frac), H, (XI, YI), method='cubic')

# Plot heatmap
plt.figure()

contour = plt.contourf(XI, YI, ZI, 50)

# Overlay actual data points
plt.scatter(std, frac, c=H, edgecolor='k')

plt.colorbar(contour, label='Shannon entropy H')

plt.xlabel('Std of inter-gene distance')
plt.ylabel('Fraction of convergent/divergent neighbours')
plt.title('Entropy plateau vs gene organisation')

plt.show()