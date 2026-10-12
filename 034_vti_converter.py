import numpy as np
import glob
import os
import vtk
from vtk.util import numpy_support


# ============================================================
# Settings
# ============================================================

NX = 32
NY = 32
NZ = 32

# Physical box lengths
# Change these to match your Fortran code
LX = 32.0
LY = 32.0
LZ = 32.0


# ============================================================
# Convert one .dat file to .vti
# ============================================================

def convert_file(filename):

    print(f"\nReading {filename}")

    # --------------------------------------------------------
    # Read Fortran output
    #
    # Columns:
    # i  j  k  phi  c  Curve  grad_u
    # --------------------------------------------------------
    data = np.loadtxt(filename)

    if data.shape != (NX * NY * NZ, 7):
        raise ValueError(
            f"{filename}: expected {(NX*NY*NZ, 7)}, "
            f"got {data.shape}"
        )

    # Integer grid coordinates
    ijk = data[:, :3].astype(int)

    # --------------------------------------------------------
    # Check coordinates
    # --------------------------------------------------------

    expected = {
        (i, j, k)
        for i in range(1, NX + 1)
        for j in range(1, NY + 1)
        for k in range(1, NZ + 1)
    }

    actual = set(map(tuple, ijk))

    missing = expected - actual

    if missing:
        raise ValueError(
            f"{filename}: missing grid points: {sorted(missing)}"
        )

    if len(actual) != len(ijk):
        raise ValueError(
            f"{filename}: duplicate grid coordinates found"
        )

    # --------------------------------------------------------
    # Extract fields
    # --------------------------------------------------------

    phi_data   = data[:, 3]
    c_data     = data[:, 4]
    curve_data = data[:, 5]
    grad_data  = data[:, 6]

    # --------------------------------------------------------
    # Physical spacing
    #
    # Your Fortran uses:
    # dx = Lx/Nx
    # dy = Ly/Ny
    # dz = Lz/Nz
    # --------------------------------------------------------

    dx = LX / NX
    dy = LY / NY
    dz = LZ / NZ

    # --------------------------------------------------------
    # Put data into [k,j,i] arrays.
    #
    # VTK uses x fastest, then y, then z.
    # --------------------------------------------------------

    phi   = np.empty((NZ, NY, NX), dtype=np.float64)
    c     = np.empty((NZ, NY, NX), dtype=np.float64)
    curve = np.empty((NZ, NY, NX), dtype=np.float64)
    grad  = np.empty((NZ, NY, NX), dtype=np.float64)

    for row in data:

        i = int(row[0]) - 1
        j = int(row[1]) - 1
        k = int(row[2]) - 1

        phi[k, j, i]   = row[3]
        c[k, j, i]     = row[4]
        curve[k, j, i] = row[5]
        grad[k, j, i]  = row[6]

    # --------------------------------------------------------
    # Create VTK ImageData
    #
    # NX*NY*NZ cells require
    # (NX+1)*(NY+1)*(NZ+1) points.
    # --------------------------------------------------------

    image = vtk.vtkImageData()

    image.SetDimensions(NX + 1, NY + 1, NZ + 1)
    image.SetSpacing(dx, dy, dz)

    # Cell centers should correspond to the original
    # Fortran grid points.
    image.SetOrigin(
        -0.5 * dx,
        -0.5 * dy,
        -0.5 * dz
    )

    # --------------------------------------------------------
    # Add cell data
    # --------------------------------------------------------

    def add_cell_array(name, array):

        vtk_array = numpy_support.numpy_to_vtk(
            array.ravel(order="C"),
            deep=True,
            array_type=vtk.VTK_DOUBLE
        )

        vtk_array.SetName(name)

        image.GetCellData().AddArray(vtk_array)

    add_cell_array("phi", phi)
    add_cell_array("c", c)
    add_cell_array("Curve", curve)
    add_cell_array("grad_u", grad)

    # Make c the active scalar
    image.GetCellData().SetActiveScalars("c")

    # --------------------------------------------------------
    # Write VTI
    # --------------------------------------------------------

    output = os.path.splitext(filename)[0] + ".vti"

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output)
    writer.SetInputData(image)

    if writer.Write() != 1:
        raise RuntimeError(f"Could not write {output}")

    print(f"  Grid       : {NX} x {NY} x {NZ}")
    print(f"  Cells      : {image.GetNumberOfCells()}")
    print(f"  Dimensions : {image.GetDimensions()}")
    print(f"  Spacing    : {dx}, {dy}, {dz}")
    print(f"  Written    : {output}")


# ============================================================
# Convert all step_*.dat files
# ============================================================

for filename in sorted(glob.glob("step_*.dat")):

    try:
        convert_file(filename)

    except Exception as e:
        print(f"ERROR converting {filename}: {e}")
