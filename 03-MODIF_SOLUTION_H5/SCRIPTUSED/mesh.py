"""
# ===================================================================================================================
#  CFD Solution Modifier - Mesh reader
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mécanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 29 June 2026
#  Last Modified   : 29 June 2026
#  Version         : 1.0.01
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  This module load the mesh.h5 regarding is structure
#
# -------------------------------------------------------------------------------------------------------------------
#  COPYRIGHT NOTICE
# -------------------------------------------------------------------------------------------------------------------
#  © 2026 Lilian CHAPUIS – All Rights Reserved.
#
#  This file and its structure are protected by intellectual property rights.
#  Unauthorized copying, distribution, modification, or use of this file
#  without prior written permission of the author is strictly prohibited.
#
# ===================================================================================================================
"""

"""
# ===================================================================================================================
#  Imports from standard library
# ===================================================================================================================
"""
import numpy as np
import h5py

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def load_mesh_coords(mesh_file):
    """
    Load nodal X, Y, Z coordinates from the HDF5 mesh file.
    Tries several common YALES2 paths automatically.
    Prints the full file structure if nothing is found.
    """
    CANDIDATE_PATHS = [
        ("Coordinates", "x", "y", "z"),
        ("Coordinates", "X", "Y", "Z"),
        ("Nodes",       "x", "y", "z"),
        ("mesh/Coordinates", "x", "y", "z"),
        ("Geometry/Nodes",   "x", "y", "z"),
        ("",            "x", "y", "z"),
        ("",            "X", "Y", "Z"),
    ]
    with h5py.File(mesh_file, "r") as f:
        for (group, dx, dy, dz) in CANDIDATE_PATHS:
            try:
                base = f[group] if group else f
                x = base[dx][:]
                y = base[dy][:]
                z = base[dz][:]
                print(f"  ✓ Coordinates found in '{group or '/'}' : {dx},{dy},{dz}")
                print(f"    {len(x):,} nodes")
                return x.astype(np.float64), y.astype(np.float64), z.astype(np.float64)
            except (KeyError, TypeError):
                continue

        print("⚠  Coordinates not found automatically. File structure:")
        def _show(name, obj):
            tag = "DS" if isinstance(obj, h5py.Dataset) else "GR"
            shp = f" {obj.shape}" if isinstance(obj, h5py.Dataset) else ""
            print(f"  [{tag}] {name}{shp}")
        f.visititems(_show)
        raise KeyError(
            "Cannot locate XYZ coordinates.\n"
            "Add the correct path to CANDIDATE_PATHS in load_mesh_coords()."
        )
