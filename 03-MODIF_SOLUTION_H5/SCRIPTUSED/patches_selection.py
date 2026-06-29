"""
# ===================================================================================================================
#  CFD Solution Modifier - Patches selection
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
#  This module store all the function to let the user select a patches dimension
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
import sys

"""
# ===================================================================================================================
#  Imports from other modules
# ===================================================================================================================
"""

from SCRIPTUSED.predefined_patches import PREDEFINED_PATCHES
from SCRIPTUSED.utils import prompt_float_opt,prompt_float,prompt_choice

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def select_patch():
    """
    Display predefined patches and manual option.
    Returns a dict with the geometric parameters of the chosen patch.
    """
    print("\n" + "═" * 55)
    print("  VOLUME DEFINITION")
    print("═" * 55)
    print("\n  Predefined patches:")

    for i, p in enumerate(PREDEFINED_PATCHES, 1):
        cyl_type = "hollow" if p["r_inner"] > 0 else "solid"
        c = p["center"]
        print(f"\n  [{i}] {p['name']}")
        print(f"       {p['description']}")
        print(f"       Axis {p['axis'].upper()}  |  center=({c[0]}, {c[1]}, {c[2]})")
        if p["r_inner"] > 0:
            print(f"       height={p['height']}  |  r_inner={p['r_inner']}  r_outer={p['r_outer']}  [{cyl_type}]")
        else:
            print(f"       height={p['height']}  |  r={p['r_outer']}  [{cyl_type}]")
        print(f"       transition={p['transition_thickness']}")

    print(f"\n  [{len(PREDEFINED_PATCHES)+1}] Manual definition")
    print()

    n_choices = len(PREDEFINED_PATCHES) + 1
    while True:
        raw = input(f"  Your choice [1–{n_choices}]: ").strip()
        try:
            choice = int(raw)
            if 1 <= choice <= n_choices:
                break
        except ValueError:
            pass
        print(f"  ⚠  Enter an integer between 1 and {n_choices}.")

    # Predefined patch
    if choice <= len(PREDEFINED_PATCHES):
        patch = PREDEFINED_PATCHES[choice - 1].copy()
        print(f"\n  ✓ Selected patch: {patch['name']}")
        edit = input("\n  Edit patch parameters? (yes/no) [no]: ").strip().lower()
        if edit in ("yes", "y"):
            patch = _edit_patch(patch)
        return patch

    # Manual definition
    return _define_manual_patch()


def _edit_patch(patch):
    """Allow the user to modify parameters of a predefined patch."""
    print("\n  Edit parameters (press Enter to keep current value):")

    cx, cy, cz = patch["center"]
    new_cx = prompt_float_opt(f"    x0 [{cx}]: ", cx)
    new_cy = prompt_float_opt(f"    y0 [{cy}]: ", cy)
    new_cz = prompt_float_opt(f"    z0 [{cz}]: ", cz)
    patch["center"] = (new_cx, new_cy, new_cz)

    raw_axis = input(f"    Axis [{patch['axis']}] (x/y/z): ").strip().lower()
    if raw_axis in ("x", "y", "z"):
        patch["axis"] = raw_axis

    patch["height"]  = prompt_float_opt(f"    Height [{patch['height']}]: ",       patch["height"])
    patch["r_outer"] = prompt_float_opt(f"    Outer radius [{patch['r_outer']}]: ", patch["r_outer"])

    cyl_type_current = "hollow" if patch["r_inner"] > 0 else "solid"
    raw_type = input(f"    Cylinder type [{cyl_type_current}] (solid/hollow): ").strip().lower()
    if raw_type == "solid":
        patch["r_inner"] = 0.0
    elif raw_type == "hollow":
        patch["r_inner"] = prompt_float_opt(
            f"    Inner radius [{patch['r_inner']}]: ", patch["r_inner"]
        )

    patch["transition_thickness"] = prompt_float_opt(
        f"    Transition thickness [{patch['transition_thickness']}]: ",
        patch["transition_thickness"]
    )
    return patch


def _define_manual_patch():
    """Full manual patch definition."""
    print("\n── Manual definition ──")
    cx = prompt_float("  x0: ")
    cy = prompt_float("  y0: ")
    cz = prompt_float("  z0: ")

    axis   = prompt_choice("  Axis (x/y/z): ", ["x", "y", "z"])
    height = prompt_float("  Height: ")

    cyl_type = prompt_choice("  Cylinder type (solid/hollow): ", ["solid", "hollow"])
    r_outer  = prompt_float("  Outer radius: ")
    r_inner  = 0.0
    if cyl_type == "hollow":
        r_inner = prompt_float("  Inner radius: ")
        if r_inner >= r_outer:
            sys.exit("✗  r_inner must be < r_outer.")

    default_thick = round(r_outer * 0.10, 6)
    thick = prompt_float_opt(
        f"  Transition thickness [default={default_thick}]: ",
        default_thick
    )

    return {
        "name":        "manual",
        "description": "Manually defined patch",
        "center":      (cx, cy, cz),
        "axis":        axis,
        "height":      height,
        "r_inner":     r_inner,
        "r_outer":     r_outer,
        "transition_thickness": thick,
    }

