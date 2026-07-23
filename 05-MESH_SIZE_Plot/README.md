# Mesh Slice Viewer

Python tool to slice a tetrahedral CFD mesh (HDF5 file, HIP/AVBP-style) along a user-defined cutting plane, and export the resulting cross-section as an image, colored by local cell size.

Developed at **IMFT** (Institut de Mécanique des Fluides de Toulouse).

---

## Table of contents

- [Features](#features)
- [Requirements](#requirements)
- [Input file format](#input-file-format)
- [Usage](#usage)
  - [Inspect an HDF5 file](#inspect-an-hdf5-file)
  - [Interactive mode](#interactive-mode)
  - [Command-line mode](#command-line-mode)
  - [Arguments table](#arguments-table)
- [Predefined zoom zones](#predefined-zoom-zones)
- [Available colormaps](#available-colormaps)
- [Output files](#output-files)
- [How the slicing works](#how-the-slicing-works)
- [Performance on large meshes](#performance-on-large-meshes)
- [Repository structure](#repository-structure)
- [Author](#author)

---

## Features

- Loads a tetrahedral mesh from an HDF5 file (separate x/y/z coordinate datasets + a flat tet→node connectivity array).
- Cutting-plane definition:
  - **axis-aligned** (`x`, `y` or `z` + position), or
  - **arbitrary plane** (point + normal vector).
- **Exact** plane/tetrahedron intersection (marching tetrahedra): each tet crossed by the plane produces a triangle or a quadrilateral, not an approximate projection. Regions with no mesh (e.g. a solid insert) therefore stay correctly empty in the figure.
- Coloring by characteristic cell size (equivalent diameter, in mm), log color scale.
- Optional zoom on predefined zones (useful for combustor geometries: swirler, flame, chamber).
- Colormap selection.
- PNG export (main plot) + PDF export (separate colorbar legend).
- Fully interactive mode (console prompts) or fully scriptable via command-line arguments.

---

## Requirements

```bash
pip install h5py numpy matplotlib
```

Python ≥ 3.8 recommended.

---

## Input file format

The script expects an HDF5 file structured as follows (HIP/AVBP-style):

```
Coordinates/
├── x        (N,)   float
├── y        (N,)   float
└── z        (N,)   float

Connectivity/
└── tet->node   (4*M,) or (M,4)   int   -- node indices, 1-based or 0-based
```

- `N` = number of nodes, `M` = number of tetrahedra.
- Connectivity can be 1-based (Fortran) or 0-based: the script detects and corrects this automatically.
- Group/dataset names can differ from the defaults, see `--points-key` and `--cells-key` below.

To check the actual structure of your file before running anything, use `--inspect` mode (see below).

---

## Usage

### Inspect an HDF5 file

Prints the full tree (groups, datasets, shapes, dtypes) without plotting anything:

```bash
python mesh_slice_viridis.py mesh.h5 --inspect
```

### Interactive mode

With no plane argument, the script asks questions step by step (mesh bounding box, axis/position or point/normal, zoom, colormap):

```bash
python mesh_slice_viridis.py mesh.h5
```

### Command-line mode

**Axis-aligned cut**, example: plane z = 0.5 m

```bash
python mesh_slice_viridis.py mesh.h5 --axis z --position 0.5
```

**Arbitrary plane**, point + normal

```bash
python mesh_slice_viridis.py mesh.h5 --point 0 0 0 --normal 1 1 0
```

**With zoom(s) and colormap forced** (no prompts asked):

```bash
python mesh_slice_viridis.py mesh.h5 --axis y --position 0.0 \
    --zone inlet flame --cmap viridis
```

Each zone listed in `--zone` generates a separate image.

**With an explicit output base name**:

```bash
python mesh_slice_viridis.py mesh.h5 --axis z --position 0.02 -o my_slice
```

**With non-default HDF5 keys**:

```bash
python mesh_slice_viridis.py mesh.h5 --axis z --position 0.5 \
    --points-key Coords --cells-key Connectivity/tetra
```

### Arguments table

| Argument | Description | Default |
|---|---|---|
| `h5_file` | Path to the input HDF5 file (positional) | — |
| `-o`, `--output` | Base name (no extension) for the output files | auto-generated |
| `--inspect` | Print the file structure and exit | `False` |
| `--points-key` | HDF5 group containing `x`/`y`/`z` | `Coordinates` |
| `--cells-key` | Path to the tet→node connectivity dataset | `Connectivity/tet->node` |
| `--axis` | `x`, `y` or `z` — axis-aligned plane | `None` |
| `--position` | Position of the plane along `--axis` | `0.0` |
| `--point` | Point on the plane (`PX PY PZ`) — arbitrary plane | `None` |
| `--normal` | Plane normal (`NX NY NZ`) — arbitrary plane | `None` |
| `--zone` | One or more zoom zones: `none inlet flame comb` | interactive prompt |
| `--cmap` | Colormap, one of the list below | interactive prompt |

> `--axis`/`--position` and `--point`/`--normal` are mutually exclusive: the script exits with an error if both methods are provided at once.

---

## Predefined zoom zones

Only meaningful for a plane aligned on `x` or `y` (where `z` is the vertical plot axis):

| CLI key | Name | Z range (m) |
|---|---|---|
| `inlet` | INLET_SWIRLER | [-0.11, 0.01] |
| `flame` | FLAME | [-0.04, 0.05] |
| `comb` | COMB_CHAMBER | [-0.01, 0.14] |
| `none` | no zoom (full range) | — |

---

## Available colormaps

`viridis` (default), `Reds_r`, `Blues_r`, `gray`, `cividis`, `BuGn_r`, `RdPu_r`, `BuPu_r`.

---

## Output files

For each requested zone:

- `<basename>.png` — the mesh cross-section colored by cell size.
- `<basename>_colorbar.pdf` — the colorbar legend alone, generated once per run (shared across all zones if several are requested).

The base name is built automatically from the chosen parameters if `-o` is not provided, e.g.:
`slice_z0p5_viridis.png`, `slice_oblique_pt0p0_0p0_n1p1_0_viridis.png`.

---

## How the slicing works

For each tetrahedron crossed by the plane, the exact intersection is computed edge by edge (6 edges per tet):

- If **1 vertex** is on one side of the plane and **3 on the other** → the plane crosses **3 edges** → the section is a **triangle**.
- If **2 vertices** are on each side → the plane crosses **4 edges** → the section is a **quadrilateral**.

The intersection points are sorted by angle around their centroid (the intersection of a plane with a convex tetrahedron is itself convex), which guarantees a correctly ordered, non self-intersecting polygon for display.

---

## Performance on large meshes

On meshes with tens to hundreds of millions of tetrahedra, this script is optimize to be as fast as possible

---

## Author

**Lilian Chapuis** — IMFT, Institut de Mécanique des Fluides de Toulouse

(c) 2026 Lilian Chapuis — All Rights Reserved. Unauthorized copying, distribution, modification, or use of this code without prior written permission of the author is strictly prohibited.
