import numpy as np
import meshio
import glob
import os

# Find all .dat files matching pattern
dat_files = sorted(glob.glob("step_*.dat"))

for file in dat_files:
    try:
        # Load data (skip header row if present)
        data = np.loadtxt(file, skiprows=1)

        # Extract coordinates (first 3 columns)
        points = data[:, :3]

        # Extract fields
        phi = data[:, 3]
        c = data[:, 4]   # second field

        # Define mesh cells (point cloud)
        cells = [("vertex", np.arange(len(points)).reshape(-1, 1))]

        # Store both fields as point data
        point_data = {
            "phi": phi,
            "c": c
        }

        # Create and write mesh
        mesh = meshio.Mesh(points=points, cells=cells, point_data=point_data)
        vtu_file = file.replace(".dat", ".vtu")
        mesh.write(vtu_file)

        print(f"Converted: {file} → {vtu_file}")

    except Exception as e:
        print(f"Error converting {file}: {e}")

print("\nConversion complete!")
