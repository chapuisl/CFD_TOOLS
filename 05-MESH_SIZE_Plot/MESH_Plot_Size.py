#!/usr/bin/env python3
"""
mesh_slice_viridis.py

Read a mesh stored in an HDF5 file (points + cell connectivity),
cut it with a user-defined plane, and save a PNG image of the
resulting cross-section, colored with the 'viridis' colormap
according to local mesh (cell) size.

------------------------------------------------------------------
HOW IT WORKS
------------------------------------------------------------------
1. The HDF5 file is opened and scanned for two datasets:
     - a "points" dataset of shape (N, 3)  -> node coordinates
     - a "cells" dataset (or group of datasets) containing the
       connectivity (indices of points forming each cell/element)

   Common dataset names are tried automatically (see
   POINTS_KEY_CANDIDATES / CELLS_KEY_CANDIDATES below). If your
   file uses different names, pass them explicitly with
   --points-key and --cells-key, or run with --inspect to print
   the full HDF5 structure first.

2. A cutting plane is defined either by:
     - an axis + a position   (--axis x|y|z --position VALUE)
     - or a point + a normal  (--point px py pz --normal nx ny nz)
   If neither is given on the command line, the script INTERACTIVELY
   ASKS the user (in the terminal) how to define the plane: axis +
   position, or point + normal, showing the mesh bounding box to
   help choose sensible values.

3. Cells whose center lies within a thin slab around the plane
   (thickness controlled by --tol) are selected. This gives the
   cross-section of the mesh.

4. For each selected cell, a characteristic size is computed
   (mean edge length of the cell). The cross-section is drawn as
   a wireframe (cell edges) projected onto the cutting plane,
   colored with the viridis colormap according to that size.

5. The figure is saved as a PNG file.

------------------------------------------------------------------
USAGE EXAMPLES
------------------------------------------------------------------
# Just inspect the HDF5 file structure first
python mesh_slice_viridis.py point.h5 --inspect

# Cut with the plane z = 0.5
python mesh_slice_viridis.py point.h5 --axis z --position 0.5 -o slice_z.png

# Cut with an arbitrary plane (point + normal)
python mesh_slice_viridis.py point.h5 --point 0 0 0 --normal 1 1 0 -o slice_oblique.png

# If dataset names are not auto-detected
python mesh_slice_viridis.py point.h5 --points-key /mesh/points --cells-key /mesh/cells
"""

import argparse
import sys

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")  
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import PolyCollection


# ----------------------------------------------------------------
# Common dataset name candidates (adjust/extend if needed)
# ----------------------------------------------------------------
POINTS_KEY_CANDIDATES = [
    "points", "Points", "coordinates", "Coordinates",
    "nodes", "Nodes", "geometry/points", "mesh/points",
    "Mesh/Points", "/Mesh/coordinates",
]

CELLS_KEY_CANDIDATES = [
    "cells", "Cells", "connectivity", "Connectivity",
    "elements", "Elements", "topology", "Topology",
    "mesh/cells", "Mesh/Cells", "/Mesh/connectivity",
]


def inspect_h5(path):
    """Print the full structure (groups, datasets, shapes) of the HDF5 file."""
    print(f"Structure of '{path}':")

    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"  Dataset: /{name}  shape={obj.shape}  dtype={obj.dtype}")
        else:
            print(f"  Group:   /{name}")

    with h5py.File(path, "r") as f:
        f.visititems(visitor)


def _find_dataset(h5file, explicit_key, candidates, expected_ndim=None):
    """Try an explicit key first, then a list of common candidate names."""
    keys_to_try = [explicit_key] if explicit_key else []
    keys_to_try += candidates

    for key in keys_to_try:
        if key is None:
            continue
        key = key.lstrip("/")
        if key in h5file:
            obj = h5file[key]
            if isinstance(obj, h5py.Dataset):
                return key, obj[()]
    return None, None


def load_mesh(path, points_key=None, cells_key=None):
    """
    Load node coordinates and cell connectivity from the HDF5 file.

    Returns
    -------
    points : (N, 3) float array
    cells  : list of 1D integer arrays (variable number of vertices per cell)
    # """


    with h5py.File(path, "r") as f:

        # ----- Lecture des coordonnées -----
        x = f["Coordinates/x"][:]
        y = f["Coordinates/y"][:]
        z = f["Coordinates/z"][:]

        points = np.column_stack((x, y, z))

        print(f"Loaded {len(points)} nodes")

        # ----- Lecture des tétraèdres -----
        tet = f["Connectivity/tet->node"][:].reshape(-1, 4)

        # Conversion éventuelle vers l'indexation Python
        if tet.min() == 1:
            tet -= 1

        cells = [c.astype(int) for c in tet]

        print(f"Loaded {len(cells)} tetrahedra")

        return points, cells


def cell_edges(cell):
    """Return list of (i, j) vertex-index pairs forming the edges of a cell."""
    n = len(cell)
    return [(cell[k], cell[(k + 1) % n]) for k in range(n)]

def tetra_volume(points, cell):
    """
    Volume of a tetrahedral cell.
    """
    a, b, c, d = points[cell[:4]]

    volume = abs(
        np.linalg.det(
            np.column_stack((b-a, c-a, d-a))
        )
    ) / 6.0

    return volume

def cell_characteristic_size(points, cell):
    # """Mean edge length of a cell -> used as the local mesh-size metric."""
    # edges = cell_edges(cell)
    # lengths = [np.linalg.norm(points[i] - points[j]) for i, j in edges]
    # return float(np.mean(lengths))
    """
    Equivalent diameter based on tetrahedron volume.
    h = diameter of sphere with same volume.
    """

    V = tetra_volume(points, cell)

    hcell = 2.0 * (3.0 * V / (4.0 * np.pi))**(1.0/3.0) *1e3

    return float(hcell)


def build_plane(axis, position, point, normal):
    """
    Build (plane_point, plane_normal, u_axis, v_axis) from user input.
    u_axis / v_axis are two orthonormal in-plane vectors, used to
    project the 3D cross-section into 2D coordinates for plotting.
    """
    if axis is not None:
        plane_point = np.zeros(3)
        plane_normal = np.zeros(3)
        idx = {"x": 0, "y": 1, "z": 2}[axis]
        plane_point[idx] = position
        plane_normal[idx] = 1.0
    else:
        plane_point = np.array(point, dtype=float)
        plane_normal = np.array(normal, dtype=float)

    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Choix de l'orientation de l'affichage
    if axis == "x":
        # Plan YZ
        u_axis = np.array([0.0, 1.0, 0.0])   # horizontal = Y
        v_axis = np.array([0.0, 0.0, 1.0])   # vertical = Z

    elif axis == "y":
        # Plan XZ
        u_axis = np.array([1.0, 0.0, 0.0])   # horizontal = X
        v_axis = np.array([0.0, 0.0, 1.0])   # vertical = Z

    elif axis == "z":
        # Plan XY
        u_axis = np.array([1.0, 0.0, 0.0])   # horizontal = X
        v_axis = np.array([0.0, 1.0, 0.0])   # vertical = Y

    else:
        # Plan arbitraire
        helper = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(helper, plane_normal)) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])

        u_axis = np.cross(plane_normal, helper)
        u_axis /= np.linalg.norm(u_axis)
        v_axis = np.cross(plane_normal, u_axis)

    return plane_point, plane_normal, u_axis, v_axis


def slice_mesh(points, cells, plane_point, plane_normal, tol,u_axis, v_axis):
    """
    Select cells whose centroid lies within `tol` of the cutting plane.

    Returns
    -------
    segments : list of ((x1, y1), (x2, y2)) edge segments in 2D plane
               coordinates, ready for a LineCollection
    values   : array of per-segment mesh-size values (for coloring)
    """

    polygons = []
    values = []

    for cell in cells:
        if len(cell) < 2:
            continue
        verts = points[cell]
        centroid = verts.mean(axis=0)
        dist = np.dot(centroid - plane_point, plane_normal)

        if abs(dist) > tol:
            continue

        size = cell_characteristic_size(points, cell)

        poly = []

        for p in verts:
            poly.append([
                np.dot(p - plane_point, u_axis),
                np.dot(p - plane_point, v_axis),
            ])

        polygons.append(poly)
        values.append(size)

    return polygons, np.array(values)


def plot_slice(polygons, values, output_path, title=None):

    if len(polygons) == 0:
        raise RuntimeError("No cells selected.")

    fig, ax = plt.subplots(figsize=(9,8))

    available_maps = {
        "1": "viridis",
        "2": "Reds_r",
        "3": "Blues_r",
        "4": "gray",
        "5": "cividis",
        "6": "bwr"
    }
    default_cmap = "Reds_r" 
    print("\nChoose a colormap for mesh size visualization:")
    for key, cmap in available_maps.items():
        print(f"{key}: {cmap}")

    max_attempts = 2
    Color_choice = False
    for attempt in range(max_attempts):
        choice = input("\nEnter your choice (1-6): ")

        if choice in available_maps:
            selected = available_maps[choice]
            print(f"Selected colormap: {selected}")
            Color_choice = True
            break

        else:
            print(f"Invalid choice. Attempt {attempt+1}/{max_attempts}")
    
    if Color_choice is False:
        selected = default_cmap

    values_log = np.maximum(values, 1e-12)

    pc = PolyCollection(
    polygons,
    array=values,
    cmap=selected,
    norm=LogNorm(
        vmin=max(values.min(),1e-12),
        vmax=values.max()
    ),
    edgecolors="black",
    linewidths=0.1,
)

    ax.add_collection(pc)

    pts = np.vstack([np.asarray(p) for p in polygons])

    # margin = 0.02 * (pts.max(axis=0) - pts.min(axis=0))

    ax.set_xlim(pts[:,0].min()-margin[0], pts[:,0].max()+margin[0])
    ax.set_ylim(pts[:,1].min()-margin[1], pts[:,1].max()+margin[1])

    ax.set_aspect("equal")

    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label("Cell size [mm]")
    cbar.ax.yaxis.set_major_formatter(
        FormatStrFormatter('%.3f')
    )

    if title:
        ax.set_title(title)

    fig.tight_layout()
    fig.savefig('Cell_size_plot', dpi=300)
    plt.close(fig)

def ask_plane_interactively(points):
    """
    Interactively ask the user how they want to define the cutting plane.

    Returns (axis, position, point, normal) in the same format expected
    by build_plane(): either (axis, position, None, None) or
    (None, None, point, normal).
    """
    bmin = points.min(axis=0)
    bmax = points.max(axis=0)
    print("\nMesh bounding box:")
    print(f"  x: [{bmin[0]:.6g}, {bmax[0]:.6g}]")
    print(f"  y: [{bmin[1]:.6g}, {bmax[1]:.6g}]")
    print(f"  z: [{bmin[2]:.6g}, {bmax[2]:.6g}]")

    print("\nHow do you want to define the cutting plane?")
    print("  1) Axis-aligned plane (choose x, y or z and a position)")
    print("  2) Arbitrary plane (a point on the plane + a normal vector)")

    while True:
        choice = input("Enter choice [1/2] (default: 1): ").strip() or "1"
        if choice in ("1", "2"):
            break
        print("Please enter 1 or 2.")

    if choice == "1":
        while True:
            axis = input("Which axis? [x/y/z]: ").strip().lower()
            if axis in ("x", "y", "z"):
                break
            print("Please enter x, y or z.")

        idx = {"x": 0, "y": 1, "z": 2}[axis]
        default_pos = float((bmin[idx] + bmax[idx]) / 2.0)
        raw = input(f"Position along {axis} (default: {default_pos:.6g}, "
                     f"range [{bmin[idx]:.6g}, {bmax[idx]:.6g}]): ").strip()
        position = float(raw) if raw else default_pos
        return axis, position, None, None

    else:
        print("Enter the coordinates of a point on the plane (default: mesh center).")
        center = (bmin + bmax) / 2.0
        raw = input(f"Point [px py pz] (default: {center[0]:.6g} {center[1]:.6g} {center[2]:.6g}): ").strip()
        if raw:
            point = [float(v) for v in raw.split()]
            if len(point) != 3:
                raise ValueError("Expected 3 values for the point (px py pz).")
        else:
            point = list(center)

        while True:
            raw = input("Normal vector [nx ny nz] (e.g. 0 0 1 for a horizontal plane): ").strip()
            normal = [float(v) for v in raw.split()] if raw else None
            if normal is not None and len(normal) == 3 and np.linalg.norm(normal) > 0:
                break
            print("Please enter 3 non-zero values, e.g. '1 0 0' or '1 1 0'.")

        return None, None, point, normal


def parse_args():
    p = argparse.ArgumentParser(
        description="Slice a mesh from an HDF5 file and export a PNG colored by mesh size."
    )
    p.add_argument("h5_file", help="Path to the input .h5 mesh file")
    p.add_argument("-o", "--output", default="mesh_slice.png", help="Output PNG path")

    p.add_argument("--inspect", action="store_true",
                    help="Print the HDF5 file structure and exit (use this to find dataset names)")

    p.add_argument("--points-key", default=None, help="Explicit HDF5 path to the points dataset")
    p.add_argument("--cells-key", default=None, help="Explicit HDF5 path to the cells dataset")

    # Plane definition: either axis+position, or point+normal
    plane_group = p.add_argument_group("cutting plane (choose ONE method)")
    plane_group.add_argument("--axis", choices=["x", "y", "z"], default=None,
                              help="Axis-aligned plane, e.g. --axis z")
    plane_group.add_argument("--position", type=float, default=0.0,
                              help="Position along --axis (default: 0.0)")
    plane_group.add_argument("--point", type=float, nargs=3, default=None,
                              metavar=("PX", "PY", "PZ"),
                              help="A point on the plane (for arbitrary plane)")
    plane_group.add_argument("--normal", type=float, nargs=3, default=None,
                              metavar=("NX", "NY", "NZ"),
                              help="Plane normal vector (for arbitrary plane)")

    p.add_argument("--tol", type=float, default=None,
                    help="Slab half-thickness around the plane used to select cells "
                         "(default: auto = 1%% of the mesh bounding-box diagonal)")

    return p.parse_args()


def main():
    args = parse_args()

    if args.inspect:
        inspect_h5(args.h5_file)
        return

    use_axis = args.axis is not None
    use_point_normal = args.point is not None and args.normal is not None
    if use_axis and use_point_normal:
        print("ERROR: specify the cutting plane with EITHER --axis/--position "
              "OR --point/--normal (not both).", file=sys.stderr)
        sys.exit(1)

    points, cells = load_mesh(args.h5_file, args.points_key, args.cells_key)

    if use_axis or use_point_normal:
        # Plane fully specified on the command line
        axis, position, point, normal = args.axis, args.position, args.point, args.normal
    else:
        # No plane given on the command line -> ask the user interactively
        axis, position, point, normal = ask_plane_interactively(points)

    plane_point, plane_normal, u_axis, v_axis = build_plane(axis, position, point, normal)

    if args.tol is None:
        bbox_diag = np.linalg.norm(points.max(axis=0) - points.min(axis=0))
        tol = 0.01 * bbox_diag
    else:
        tol = args.tol
    print(f"Using slab half-thickness (tol) = {tol:.6g}")

    segments, values = slice_mesh(points, cells, plane_point, plane_normal, tol,u_axis, v_axis)
    print(f"Selected {len(segments)} edges for the cross-section")

    if axis is not None:
        title = f"Mesh cross-section at {axis} = {position}"
    else:
        title = f"Mesh cross-section (point={point}, normal={normal})"

    plot_slice(segments, values, args.output, title=title)


if __name__ == "__main__":
    main()