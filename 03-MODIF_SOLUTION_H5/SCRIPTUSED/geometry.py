"""
# ===================================================================================================================
#  CFD Solution Modifier - Geometry Utilities
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
#  This module provides all geometric utilities required to define cylindrical
#  regions inside the computational domain.
#
#  The implemented functions allow:
#      - Radial distance computation.
#      - Selection of nodes located inside a cylindrical volume.
#      - Construction of a smooth transition region (implemented below).
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
#  Import library
# ===================================================================================================================
"""
import numpy as np

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def radial_distance(x, y, z, center, axis):
    """Radial distance from the cylinder axis and coordinate along axis."""
    cx, cy, cz = center
    if axis == "z":
        return np.sqrt((x - cx)**2 + (y - cy)**2), z - cz
    elif axis == "x":
        return np.sqrt((y - cy)**2 + (z - cz)**2), x - cx
    elif axis == "y":
        return np.sqrt((x - cx)**2 + (z - cz)**2), y - cy
    raise ValueError(f"Invalid axis: '{axis}'. Choose x, y or z.")


def select_nodes_cylinder(x, y, z, center, axis, r_inner, r_outer, height, h_offset=0.0):
    """
    Return indices of nodes inside the defined cylinder.

    Parameters
    ----------
    center   : (x0, y0, z0) — base centre of the cylinder
    axis     : 'x', 'y' or 'z'
    r_inner  : inner radius (0 for solid cylinder)
    r_outer  : outer radius
    height   : cylinder length along the axis
    h_offset : offset of the base along the axis (default 0)
    """
    r, along = radial_distance(x, y, z, center, axis)
    mask = (
        (r     >= r_inner) &
        (r     <= r_outer) &
        (along >= h_offset) &
        (along <= h_offset + height)
    )
    return np.where(mask)[0]



def build_transition_zone(x, y, z, center, axis, r_inner, r_outer, height,
                          h_offset=0.0, transition_thickness=None):
    """
    Compute the transition zone around the cylinder boundary.

    For each node outside the core but within the transition envelope:

        alpha = distance_to_core_boundary / transition_thickness

        alpha = 0  →  core boundary  →  imposed value
        alpha = 1  →  outer limit    →  original value preserved

    The law is LINEAR and CONTINUOUS, computed per node from the
    actual distance to the boundary — no layer discretisation,
    no staircase artefacts.

    Hollow cylinder fix:
      - The inner hole (r < r_inner) is explicitly excluded from
        the transition zone (was the source of the hollow-cylinder bug).
      - A transition zone is also built around the inner boundary.

    Parameters
    ----------
    transition_thickness : total thickness of the transition zone.
                           Default = r_outer * 0.10

    Returns
    -------
    idx_trans : np.ndarray  — indices of transition nodes
    alpha     : np.ndarray  — per-node linear weight in [0, 1]
    """
    if transition_thickness is None:
        transition_thickness = r_outer * 0.10

    th = transition_thickness
    r, along = radial_distance(x, y, z, center, axis)
    h_lo = h_offset
    h_hi = h_offset + height

    # Core mask (nodes already handled by modify_species)
    in_core = (
        (r     >= r_inner) & (r     <= r_outer) &
        (along >= h_lo)    & (along <= h_hi)
    )

    # Inner hole mask (hollow cylinder) — must NOT be modified
    in_hole = (r < r_inner) & (along >= h_lo - th) & (along <= h_hi + th)

    # Transition zone around the OUTER boundary
    in_trans = (
        (r     <= r_outer + th) &
        (along >= h_lo - th) & (along <= h_hi + th) &
        ~in_core & ~in_hole
    )

    # Transition zone around the INNER boundary (hollow cylinder only)
    if r_inner > 0.0:
        in_inner_trans = (
            (r >= max(0.0, r_inner - th)) & (r < r_inner) &
            (along >= h_lo - th) & (along <= h_hi + th) &
            ~in_hole
        )
        in_trans = in_trans | in_inner_trans

    idx_trans = np.where(in_trans)[0]
    if len(idx_trans) == 0:
        return idx_trans, np.array([], dtype=np.float64)

    r_t     = r[idx_trans]
    along_t = along[idx_trans]

    # Signed distance to each core boundary (positive = outside core)
    d_r_ext = r_t - r_outer                    # outer radial boundary
    d_h_bot = h_lo - along_t                   # axial bottom boundary
    d_h_top = along_t - h_hi                   # axial top boundary
    d_r_int = (r_inner - r_t) if r_inner > 0.0 else np.full_like(r_t, -np.inf)

    # Distance to core = max of signed distances (≥ 0 since nodes are outside)
    d_to_core = np.maximum.reduce([d_r_ext, d_h_bot, d_h_top, d_r_int])
    d_to_core = np.clip(d_to_core, 0.0, th)

    # Linear law: 0 at core boundary → 1 at outer transition limit
    alpha = d_to_core / th

    return idx_trans, alpha
