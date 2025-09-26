import numpy as np
import matplotlib.pyplot as plt
import imageio
import glob
import os

# Collect all files step_*.dat
files = sorted(glob.glob("step_*.dat"), key=lambda x: int(x.split("_")[1].split(".")[0]))

# First pass: compute global axis and phi limits
xmin = ymin = zmin = float("inf")
xmax = ymax = zmax = float("-inf")
phi_min = float("inf")
phi_max = float("-inf")

for f in files:
    data = np.loadtxt(f)
    x, y, z, grad_phi, phi = data[:,0], data[:,1], data[:,2], data[:,4], data[:,5]
    mask = np.abs(phi) > 0.7
    if np.any(mask):
        x, y, z, grad_phi, phi = x[mask], y[mask], z[mask], grad_phi[mask], phi[mask]
        xmin, ymin, zmin = 0,0,0
        xmax, ymax, zmax = 16,16,16
        phi_min = min(phi_min, phi.min())
        phi_max = max(phi_max, phi.max())

# Second pass: plot each frame with fixed axis + fixed color scale
frames = []
for f in files:
    data = np.loadtxt(f)
    x, y, z, phi = data[:,0], data[:,1], data[:,2], data[:,5]
    mask = np.abs(phi) > 0.7
    x, y, z, phi = x[mask], y[mask], z[mask], phi[mask]
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(
        x, y, z, c=phi,
        cmap="viridis",
        s=40, marker="o",
        vmin=phi_min, vmax=phi_max   # fixed color scale
    )
    
    # Fixed axis scaling (global)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_zlim([zmin, zmax])
    
    ax.set_xlabel("Grid X")
    ax.set_ylabel("Grid Y")
    ax.set_zlabel("Grid Z")
    ax.set_title(f"Step {f.split('_')[1].split('.')[0]}")
    
    fig.colorbar(p, ax=ax, label="Phi")
    
    # Save temp frame
    fname = f"frame_{f.split('_')[1].split('.')[0]}.png"
    plt.savefig(fname, dpi=100)
    plt.close(fig)
    
    frames.append(imageio.imread(fname))

# Save animation
imageio.mimsave("movie.gif", frames, fps=10)   # GIF
imageio.mimsave("movie.mp4", frames, fps=10)   # MP4

# Cleanup PNGs if you want
for f in glob.glob("frame_*.png"):
    os.remove(f)

