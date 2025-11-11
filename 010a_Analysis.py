import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

# ---------------------------------------------------------
# USER INPUT
# ---------------------------------------------------------
Nx, Ny, Nz = 16, 16, 16          # Must match Fortran
directory = "."                  # Folder containing step_*.dat
output_video = "hernandez_3d.mp4"
dt = 0.01                        # Must match Fortran time increment
# ---------------------------------------------------------

# ---------------------------------------------------------
# Collect only valid step_<integer>.dat files
# ---------------------------------------------------------
files = [
    f for f in glob.glob(os.path.join(directory, "step_*.dat"))
    if f.split("_")[1].split(".")[0].isdigit()
]

files = sorted(files, key=lambda x: int(x.split("_")[1].split(".")[0]))

# ---------------------------------------------------------
# Center slice indices
# ---------------------------------------------------------
x_mid = Ny // 2
z_mid = Nz // 2

# ---------------------------------------------------------
# Scan all data to determine fixed axis range
# ---------------------------------------------------------
c_min, c_max = 1e9, -1e9
for f in files:
    data = np.loadtxt(f).reshape(Nx, Ny, Nz, 5)
    c_line = data[x_mid, :, z_mid, 4]
    c_min = min(c_min, np.min(c_line))
    c_max = max(c_max, np.max(c_line))

# ---------------------------------------------------------
# Setup figure
# ---------------------------------------------------------
fig, ax = plt.subplots()
line, = ax.plot([], [], marker="o")

ax.set_xlim(1, Nx)
ax.set_ylim(c_min, c_max)
ax.set_xlabel("x index")
ax.set_ylabel("Concentration c")

# ---------------------------------------------------------
# Init function (prevents duplicate step_0.dat execution)
# ---------------------------------------------------------
def init():
    line.set_data([], [])
    return line,

# ---------------------------------------------------------
# Update each animation frame
# ---------------------------------------------------------
def update(frame):
    datafile = files[frame]
#    print(f"Processing: {datafile}")

    data = np.loadtxt(datafile).reshape(Nx, Ny, Nz, 5)
    c_line = data[x_mid, :, z_mid, 4]

    line.set_data(range(1, Nx + 1), c_line)

    step_num = int(os.path.basename(datafile).split("_")[1].split(".")[0])


    ax.set_title(f"Time = {step_num}")

    return line,

# ---------------------------------------------------------
# Create animation (blit + init fixes double-call issue)
# ---------------------------------------------------------
ani = FuncAnimation(
    fig,
    update,
    frames=len(files),
    init_func=init,
    interval=200,
    blit=True
)

ani.save(output_video, dpi=150, fps=5)

sum_values = []

for f in files:
    data = np.loadtxt(f).reshape(Nx, Ny, Nz, 5)
    
    # 5th column = index 4 (0-based indexing)
    c_sum = np.sum(data[:, :, :, 4])
    step_num = int(os.path.basename(f).split("_")[1].split(".")[0])

    sum_values.append((step_num, c_sum))

    print(f"Step {step_num:6d}:  Sum of c = {c_sum:.6f}")
