import numpy as np
import matplotlib.pyplot as plt

# Load a snapshot
data = np.loadtxt("step_0.dat")  # or any step

# Extract columns
x, y, z, phi = data[:,0], data[:,1], data[:,2], data[:,4]

# Pick slice y=0, z=0
mask = (y == 9) & (z == 9)  # assuming y=0, z=0 are min values
x_slice = x[mask]
phi_slice = phi[mask]

# Sort by x for plotting
sort_idx = np.argsort(x_slice)
x_slice = x_slice[sort_idx]
phi_slice = phi_slice[sort_idx]

# Plot
plt.figure(figsize=(8,5))
plt.plot(x_slice, phi_slice, '-o')
plt.xlabel("x")
plt.ylabel("phi")
plt.title("Phi along y=0, z=0")
plt.grid(True)
plt.show()
