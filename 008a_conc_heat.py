import numpy as np
import matplotlib.pyplot as plt
import imageio
import glob
import os

# Collect all files step_*.dat
files = sorted(glob.glob("step_*.dat"), key=lambda x: int(x.split("_")[1].split(".")[0]))

# First pass: compute global axis and conc limits
xmin = ymin = zmin = float("inf")
xmax = ymax = zmax = float("-inf")
conc_min = float("inf")
conc_max = float("-inf")

for f in files:
    data = np.loadtxt(f)
    x, y, z, phi, conc = data[:,0], data[:,1], data[:,2], data[:,3],data[:,4]
    mask = np.abs(conc)<1-1e-3
    if np.any(mask):
        x, y, z, conc = x[mask], y[mask], z[mask], conc[mask]
        xmin, ymin, zmin = 0,0,0
        xmax, ymax, zmax = 16,16,16
        conc_min = min(conc_min, conc.min())
        conc_max = max(conc_max, conc.max())

# Second pass: plot each frame with fixed axis + fixed color scale
frames = []
for f in files:
    data = np.loadtxt(f)
    x, y, z, phi, conc = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    mask = np.abs(conc) < 1-1e-3
    x, y, z, conc = x[mask], y[mask], z[mask], conc[mask]
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(
        x, y, z, c=conc,
        cmap="viridis",
        s=40, marker="o",
        vmin=conc_min, vmax=conc_max   # fixed color scale
    )
    
    # Fixed axis scaling (global)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_zlim([zmin, zmax])
    
    ax.set_xlabel("Grid X")
    ax.set_ylabel("Grid Y")
    ax.set_zlabel("Grid Z")
    ax.set_title(f"Step {f.split('_')[1].split('.')[0]}")
    
    fig.colorbar(p, ax=ax, label="conc")
    
    # Save temp frame
    fname = f"frame_{f.split('_')[1].split('.')[0]}.png"
    plt.savefig(fname, dpi=100)
    plt.close(fig)
    
    frames.append(imageio.imread(fname))

# Save animation
imageio.mimsave("08a_conc_movie.gif", frames, fps=10)   # GIF
imageio.mimsave("08a_conc_movie.mp4", frames, fps=10)   # MP4

# Cleanup PNGs if you want
for f in glob.glob("frame_*.png"):
    os.remove(f)

