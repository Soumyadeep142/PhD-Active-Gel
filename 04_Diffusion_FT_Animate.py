import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob

# === Load all files matching u_tXXXX.dat ===
file_list = sorted(glob.glob("u_t*.dat"))

# Read all data into a list
snapshots = []
x = None
for fname in file_list:
    data = np.loadtxt(fname)
    if x is None:
        x = data[:, 0]
    snapshots.append(data[:, 1])

snapshots = np.array(snapshots)  # shape: (nframes, N)
times = [int(f.split("_t")[1].split(".dat")[0]) for f in file_list]

# === Plot setup ===
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.set_xlim(x.min(), x.max())
ax.set_ylim(snapshots.min()*1.1, snapshots.max()*1.1)
ax.set_xlabel("x")
ax.set_ylabel("u(x, t)")

# === Animation update function ===
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def update(frame):
    line.set_data(x, snapshots[frame])
    time_text.set_text(f"t = {times[frame]}")
    return line, time_text

ani = FuncAnimation(fig, update, frames=len(snapshots),
                    init_func=init, blit=True, interval=200)

# === Save to MP4 ===
ani.save("heat_diffusion.mp4", fps=5, dpi=150)

plt.show()

