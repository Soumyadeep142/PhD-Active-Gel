import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Parameters ---
N = 128                  # Grid points
L = 2.0 * np.pi           # Domain length
kappa = 0.1               # Diffusion coefficient
NT = 500                  # Number of time steps
dt = 0.01                 # Time step size
snapint = 10              # Snapshot interval

# --- Derived ---
x = np.linspace(0, L, N, endpoint=False)
dx = L / N

# Wavenumbers (FFT ordering)
kvec = np.fft.fftfreq(N, d=dx) * 2.0 * np.pi

# --- Initial condition: Gaussian concentration at center ---
u = np.exp(- (x - L/2.0)**2 / 0.1)

# --- Time marching ---
snapshots = [u.copy()]
times = [0.0]

uh = np.fft.fft(u)  # Fourier coefficients

for n in range(1, NT+1):
    # Exact Fourier-space update for diffusion: e^{-kappa * k^2 * dt}
    uh *= np.exp(-kappa * kvec**2 * dt)
    
    # Back to physical space
    u = np.real(np.fft.ifft(uh))
    
    # Store snapshot
    if n % snapint == 0:
        snapshots.append(u.copy())
        times.append(n * dt)

# --- Convert snapshots to array ---
U_snap = np.array(snapshots).T  # Shape: (space, time_snap)

# --- Animation ---
fig, ax = plt.subplots()
line, = ax.plot(x, U_snap[:, 0], lw=2)
ax.set_ylim(0, U_snap.max()*1.1)
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")
title = ax.set_title(f"t = {times[0]:.2f}")

def update(frame):
    line.set_ydata(U_snap[:, frame])
    title.set_text(f"t = {times[frame]:.2f}")
    return line, title

ani = animation.FuncAnimation(fig, update, frames=len(times), blit=False, interval=100)
ani.save("diffusion_fft_python.mp4", fps=10, dpi=200)

plt.show()

