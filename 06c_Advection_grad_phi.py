import numpy as np
import matplotlib.pyplot as plt

# Load data: columns = x y z col4(col=grad_phi) col5 col6(phi)
data = np.loadtxt("step_0.dat")
x_all, y_all, z_all, grad_phi_all = data[:,0], data[:,1], data[:,2], data[:,6]

# Select points along y=0, z=0
mask = (y_all == 9) & (z_all == 9)
x_line = x_all[mask]
grad_phi_line = grad_phi_all[mask]

# Sort by x
sort_idx = np.argsort(x_line)
x_line = x_line[sort_idx]
grad_phi_line = grad_phi_line[sort_idx]

# Print or plot
for xi, gi in zip(x_line, grad_phi_line):
    print(f"x = {xi:.3f}, grad_phi = {gi:.6f}")

# Optional: plot
print(x_line)
print(grad_phi_line)
plt.plot(x_line, grad_phi_line, marker='o')
plt.xlabel('x')
plt.ylabel('grad_phi')
plt.title('grad_phi along y=0, z=0')
plt.show()

