"""
# ===================================================================================================================
#  Mesh Slice Viewer - Main Program
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mecanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 22 July 2026
#  Last Modified   : 22 July 2026
#  Version         : 1.0.00
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  Slice a CFD tetrahedral mesh stored in an HDF5 file (HIP/AVBP-style: separate
#  Coordinates/x,y,z datasets + a flat Connectivity/tet->node array) along a
#  user-defined plane, and export the resulting cross-section as a PNG image,
#  colored by local cell size (equivalent tetrahedron diameter).
#
#  This script orchestrates the complete workflow:
#      - HDF5 mesh inspection and loading
#      - Interactive or command-line cutting-plane definition (axis-aligned
#        or arbitrary point + normal)
#      - Optional Z-zoom on predefined combustor zones (INLET_SWIRLER,
#        FLAME, COMB_CHAMBER)
#      - Cell selection near the cutting plane and characteristic-size
#        computation
#      - Colormap selection and cross-section plotting (log color scale,
#        plain-number tick labels)
#      - PNG export of the cross-section + PDF export of the colorbar legend
#
# -------------------------------------------------------------------------------------------------------------------
#  USAGE EXAMPLES
# -------------------------------------------------------------------------------------------------------------------
#  python mesh_slice_viridis.py mesh.h5 --inspect
#  python mesh_slice_viridis.py mesh.h5 --axis z --position 0.5
#  python mesh_slice_viridis.py mesh.h5 --point 0 0 0 --normal 1 1 0
#  python mesh_slice_viridis.py mesh.h5   # fully interactive
#
# -------------------------------------------------------------------------------------------------------------------
#  COPYRIGHT NOTICE
# -------------------------------------------------------------------------------------------------------------------
#  (c) 2026 Lilian CHAPUIS - All Rights Reserved.
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
import argparse
import sys

"""
# ===================================================================================================================
#  Imports from third-party libraries
# ===================================================================================================================
"""
import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, FuncFormatter, NullLocator, NullFormatter
from matplotlib.collections import PolyCollection


"""
# ===================================================================================================================
#  Configuration / constants
# ===================================================================================================================
"""

# Predefined Z zoom zones (only meaningful for x- or y-aligned planes)
ZONES = {
    "1": ("INLET_SWIRLER", -0.115, 0.005),
    "2": ("FLAME", -0.04, 0.05),
    "3": ("COMB_CHAMBER", -0.01, 0.14),
}

# Available colormaps to choose from interactively
AVAILABLE_CMAPS = {
    "1": "viridis",
    "2": "Reds_r",
    "3": "Blues_r",
    "4": "gray",
    "5": "cividis",
    "6": "bwr",
}
DEFAULT_CMAP = "viridis"

# Maps the --zone CLI argument to the corresponding ZONES key
ZONE_CLI_MAP = {"inlet": "1", "flame": "2", "comb": "3", "none": None}


"""
# ===================================================================================================================
#  HDF5 inspection & mesh loading
# ===================================================================================================================
"""


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


def load_mesh(path, points_key=None, cells_key=None):
    """
    Load node coordinates and tetrahedral connectivity from the HDF5 file.

    Returns
    -------
    points : (N, 3) float array
    cells  : (M, 4) int array of node indices (0-based), one row per tet
    """
    points_group = points_key or "Coordinates"
    px = f"{points_group}/x"
    py = f"{points_group}/y"
    pz = f"{points_group}/z"
    ck = cells_key or "Connectivity/tet->node"

    with h5py.File(path, "r") as f:
        x = f[px][:]
        y = f[py][:]
        z = f[pz][:]
        points = np.column_stack((x, y, z))
        print(f"Loaded {len(points)} nodes from '{px}', '{py}', '{pz}'")

        tet = f[ck][:].reshape(-1, 4).astype(np.int64)
        if tet.min() == 1:
            tet -= 1  # convert from 1-based (Fortran/HIP) to 0-based indexing
        print(f"Loaded {len(tet)} tetrahedra from '{ck}'")

        return points, tet


"""
# ===================================================================================================================
#  Cell geometry utilities
# ===================================================================================================================
"""


def tetra_volume(points, cell):
    """Volume of a tetrahedral cell."""
    a, b, c, d = points[cell[:4]]
    return abs(np.linalg.det(np.column_stack((b - a, c - a, d - a)))) / 6.0


def cell_characteristic_size(points, cell):
    """
    Equivalent diameter based on tetrahedron volume.
    h = diameter of the sphere with the same volume, in mm
    (assumes mesh coordinates are in meters).
    """
    V = tetra_volume(points, cell)
    return float(2.0 * (3.0 * V / (4.0 * np.pi)) ** (1.0 / 3.0) * 1e3)


"""
# ===================================================================================================================
#  Cutting-plane geometry
# ===================================================================================================================
"""


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

    if axis == "x":
        u_axis = np.array([0.0, 1.0, 0.0])  # horizontal = Y
        v_axis = np.array([0.0, 0.0, 1.0])  # vertical   = Z
    elif axis == "y":
        u_axis = np.array([1.0, 0.0, 0.0])  # horizontal = X
        v_axis = np.array([0.0, 0.0, 1.0])  # vertical   = Z
    elif axis == "z":
        u_axis = np.array([1.0, 0.0, 0.0])  # horizontal = X
        v_axis = np.array([0.0, 1.0, 0.0])  # vertical   = Y
    else:
        helper = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(helper, plane_normal)) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        u_axis = np.cross(plane_normal, helper)
        u_axis /= np.linalg.norm(u_axis)
        v_axis = np.cross(plane_normal, u_axis)

    return plane_point, plane_normal, u_axis, v_axis


def slice_mesh(points, cells, plane_point, plane_normal, tol, u_axis, v_axis):
    """
    Select cells whose centroid lies within `tol` of the cutting plane.
    Vectorized over the (M, 4) cells array for performance on large meshes.

    Returns
    -------
    polygons : list of (4, 2) arrays -> projected cell vertices, ready
               for a PolyCollection
    values   : array of per-cell characteristic-size values (for coloring)
    """
    cell_pts = points[cells]                      # (M, 4, 3)
    centroids = cell_pts.mean(axis=1)              # (M, 3)
    dist = (centroids - plane_point) @ plane_normal
    mask = np.abs(dist) <= tol

    selected = cells[mask]
    print(f"Selected {len(selected)} cells near the cutting plane")

    p0u = plane_point @ u_axis
    p0v = plane_point @ v_axis

    polygons = []
    values = np.empty(len(selected))
    for n, cell in enumerate(selected):
        verts = points[cell]
        uv = np.column_stack((verts @ u_axis - p0u, verts @ v_axis - p0v))
        polygons.append(uv)
        values[n] = cell_characteristic_size(points, cell)

    return polygons, values


"""
# ===================================================================================================================
#  Interactive prompts
# ===================================================================================================================
"""


def ask_plane_interactively(points):
    """Ask the user how to define the cutting plane."""
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


def ask_zoom_interactively(axis):
    """
    Ask the user whether to zoom on one or several predefined Z zones.
    Only meaningful for x- or y-aligned planes, where Z is the vertical
    plot axis. One image will be generated per selected zone.

    Returns a list of zoom zones, where each entry is either a
    (name, zmin, zmax) tuple, or None for "no zoom" (full range).
    """
    if axis not in ("x", "y"):
        return [None]

    print("\nZoom on specific Z zone(s)? (only relevant for x/y cutting planes)")
    print("  0) No zoom (full Z range)")
    for key, (name, zmin, zmax) in ZONES.items():
        print(f"  {key}) {name}  (z in [{zmin}, {zmax}])")
    print("You can select several at once, e.g. '0 1 2' or '1 3', separated by spaces.")
    print("One output image will be generated per selection.")

    while True:
        raw = input("Enter choice(s) [0/1/2/3] (default: 0): ").strip() or "0"
        tokens = raw.split()
        if tokens and all(t in ("0", "1", "2", "3") for t in tokens):
            break
        print("Please enter one or more of 0, 1, 2, 3 separated by spaces.")

    zones = []
    seen = set()
    for t in tokens:
        if t in seen:
            continue
        seen.add(t)
        zones.append(None if t == "0" else ZONES[t])

    return zones


def ask_colormap_interactively():
    """Ask the user which colormap to use. Falls back to DEFAULT_CMAP."""
    print("\nChoose a colormap for the cell-size visualization:")
    for key, cmap in AVAILABLE_CMAPS.items():
        print(f"  {key}) {cmap}")

    max_attempts = 2
    for attempt in range(max_attempts):
        choice = input(f"Enter your choice (1-{len(AVAILABLE_CMAPS)}) "
                       f"(default: {DEFAULT_CMAP}): ").strip()
        if not choice:
            return DEFAULT_CMAP
        if choice in AVAILABLE_CMAPS:
            return AVAILABLE_CMAPS[choice]
        print(f"Invalid choice. Attempt {attempt + 1}/{max_attempts}")

    return DEFAULT_CMAP


"""
# ===================================================================================================================
#  Output naming
# ===================================================================================================================
"""


def sanitize(value):
    """Turn a number into a filename-safe token, e.g. -0.12 -> 'm0p12'."""
    s = f"{value:.4g}"
    return s.replace("-", "m").replace(".", "p")


def build_basename(axis, position, point, normal, zone, cmap_name):
    """Build an output file base name (no extension) from user choices."""
    if axis is not None:
        base = f"slice_{axis}{sanitize(position)}"
    else:
        base = (f"slice_oblique_pt{sanitize(point[0])}_{sanitize(point[1])}_{sanitize(point[2])}"
                f"_n{sanitize(normal[0])}_{sanitize(normal[1])}_{sanitize(normal[2])}")

    if zone is not None:
        base += f"_{zone[0]}"

    base += f"_{cmap_name}"
    return base


"""
# ===================================================================================================================
#  Plotting
# ===================================================================================================================
"""


def plot_slice(polygons, values, cmap_name, output_png, output_colorbar_pdf,
               title=None, zoom_zone=None, tick_fontsize=16):
    """Draw the cross-section as filled polygons and save PNG + colorbar PDF."""
    if len(polygons) == 0:
        raise RuntimeError(
            "No cells found near the cutting plane. "
            "Try increasing --tol or check the plane definition."
        )

    values = np.asarray(values)
    values_safe = np.maximum(values, 1e-12)

    norm = LogNorm(vmin=values_safe.min(), vmax=values_safe.max())

    fig, ax = plt.subplots(figsize=(9, 8))

    pc = PolyCollection(
        polygons,
        array=values,
        cmap=cmap_name,
        norm=norm,
        edgecolors="black",
        linewidths=0.1,
    )
    ax.add_collection(pc)

    pts = np.vstack(polygons)
    margin = 0.02 * (pts.max(axis=0) - pts.min(axis=0) + 1e-12)
    ax.set_xlim(pts[:, 0].min() - margin[0], pts[:, 0].max() + margin[0])

    if zoom_zone is not None:
        _, zmin, zmax = zoom_zone
        zmargin = 0.02 * (zmax - zmin)
        ax.set_ylim(zmin - zmargin, zmax + zmargin)
    else:
        ax.set_ylim(pts[:, 1].min() - margin[1], pts[:, 1].max() + margin[1])

    ax.set_aspect("equal")
    ax.tick_params(axis="both", which="major", labelsize=tick_fontsize, length=7, width=1.2)

    if title:
        ax.set_title(title, fontsize=tick_fontsize + 2)

    # ---- Colorbar on the main figure: plain (non power-of-10) tick labels ----
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label("Cell size [mm]", fontsize=tick_fontsize)
    cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
    cbar.ax.yaxis.set_minor_locator(NullLocator())
    cbar.ax.yaxis.set_minor_formatter(NullFormatter())
    cbar.ax.tick_params(labelsize=tick_fontsize)

    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)
    print(f"Saved cross-section image to '{output_png}'")

    # ---- Separate figure containing only the colorbar, saved as PDF ----
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap_name)
    mappable.set_array(values)

    fig_cb = plt.figure(figsize=(2.2, 8))
    ax_cb = fig_cb.add_axes([0.35, 0.05, 0.25, 0.9])
    cbar2 = fig_cb.colorbar(mappable, cax=ax_cb)
    cbar2.set_label("Cell size [mm]", fontsize=tick_fontsize)
    cbar2.ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
    cbar2.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3g}"))
    cbar2.ax.yaxis.set_minor_locator(NullLocator())
    cbar2.ax.yaxis.set_minor_formatter(NullFormatter())
    cbar2.ax.tick_params(labelsize=tick_fontsize)
    if output_colorbar_pdf is not None:
        fig_cb.savefig(output_colorbar_pdf)
        plt.close(fig_cb)
        print(f"Saved colorbar legend to '{output_colorbar_pdf}'")


"""
# ===================================================================================================================
#  Command-line interface
# ===================================================================================================================
"""


def parse_args():
    p = argparse.ArgumentParser(
        description="Slice a mesh from an HDF5 file and export a PNG colored by cell size."
    )
    p.add_argument("h5_file", help="Path to the input .h5 mesh file")
    p.add_argument("-o", "--output", default=None,
                   help="Base name (no extension) for the output PNG/PDF. "
                        "If omitted, a name is built automatically from your choices.")

    p.add_argument("--inspect", action="store_true",
                   help="Print the HDF5 file structure and exit")

    p.add_argument("--points-key", default=None,
                   help="HDF5 group containing x/y/z coordinate datasets (default: 'Coordinates')")
    p.add_argument("--cells-key", default=None,
                   help="HDF5 path to the flat tet connectivity dataset "
                        "(default: 'Connectivity/tet->node')")

    plane_group = p.add_argument_group("cutting plane (choose ONE method)")
    plane_group.add_argument("--axis", choices=["x", "y", "z"], default=None)
    plane_group.add_argument("--position", type=float, default=0.0)
    plane_group.add_argument("--point", type=float, nargs=3, default=None, metavar=("PX", "PY", "PZ"))
    plane_group.add_argument("--normal", type=float, nargs=3, default=None, metavar=("NX", "NY", "NZ"))

    p.add_argument("--zone", choices=["none", "inlet", "flame", "comb"], default=None, nargs="+",
                   help="One or more Z zoom zones for x/y-aligned planes, e.g. "
                        "'--zone inlet flame'. Skips the interactive prompt. "
                        "One output image is generated per zone.")
    p.add_argument("--cmap", default=None, choices=list(AVAILABLE_CMAPS.values()),
                   help="Colormap to use (skips the interactive prompt)")

    p.add_argument("--tol", type=float, default=None,
                   help="Slab half-thickness around the plane (default: 1%% of bbox diagonal)")

    return p.parse_args()


"""
# ===================================================================================================================
#  MAIN Function
# ===================================================================================================================
"""


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
        axis, position, point, normal = args.axis, args.position, args.point, args.normal
    else:
        axis, position, point, normal = ask_plane_interactively(points)

    plane_point, plane_normal, u_axis, v_axis = build_plane(axis, position, point, normal)

    # ---- Colormap (asked first) ----
    cmap_name = args.cmap or ask_colormap_interactively()

    # ---- Optional Z zoom(s) (only offered for x/y-aligned planes) ----
    if args.zone is not None:
        zones = []
        for z in dict.fromkeys(args.zone):  # de-duplicate, keep order
            zones.append(None if z == "none" else ZONES[ZONE_CLI_MAP[z]])
    else:
        zones = ask_zoom_interactively(axis)

    if args.tol is None:
        bbox_diag = np.linalg.norm(points.max(axis=0) - points.min(axis=0))
        tol = 0.01 * bbox_diag
    else:
        tol = args.tol
    print(f"Using slab half-thickness (tol) = {tol:.6g}")

    polygons, values = slice_mesh(points, cells, plane_point, plane_normal, tol, u_axis, v_axis)

    if axis is not None:
        base_title = f"Mesh cross-section at {axis} "
    else:
        base_title = f"Mesh cross-section (point={point}, normal={normal})"
    Colormapspdf = True
    for zone in zones:
        
        title = base_title if zone is None else f"{base_title}  [{zone[0]}]"

        basename = args.output or build_basename(axis, position, point, normal, zone, cmap_name)
        if args.output and len(zones) > 1:
            # avoid overwriting the same file when several zones are requested
            # together with an explicit --output base name
            basename += "_full" if zone is None else f"_{zone[0]}"
        output_png = f"{basename}.png"
        if Colormapspdf is True:
            output_colorbar_pdf = f"{basename}_colorbar.pdf"
        else:
            output_colorbar_pdf = None

        plot_slice(polygons, values, cmap_name, output_png, output_colorbar_pdf,
                   title=title, zoom_zone=zone)
        Colormapspdf = False


if __name__ == "__main__":
    main()