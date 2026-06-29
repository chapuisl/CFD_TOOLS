"""
# ===================================================================================================================
#  CFD Solution Modifier - Main Program
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mécanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 10 January 2026
#  Last Modified   : 17 February 2026
#  Version         : 1.0.02
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  Main entry point of the CFD solution modifier.
#
#  This script orchestrates the complete workflow:
#      - Input/output file selection
#      - Mesh loading
#      - Patch definition
#      - Node selection
#      - Species modification
#      - Solution verification
#      - XMF export
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
import shutil
import h5py
import sys
import os

"""
# ===================================================================================================================
#  Imports from other modules
# ===================================================================================================================
"""

from SCRIPTUSED.mesh import load_mesh_coords
from SCRIPTUSED.geometry import (select_nodes_cylinder, build_transition_zone)
from SCRIPTUSED.patches_selection import select_patch
from SCRIPTUSED.utils import prompt_float,prompt_choice,verify_solution
from SCRIPTUSED.xmf_writer import generate_xmf
from SCRIPTUSED.species import modify_species

"""
# ===================================================================================================================
#  MAIN Function
# ===================================================================================================================
"""

def main():
    print("=" * 60)
    print("   CFD Solution Modifier — Species & Volume")
    print("=" * 60)

    # ── Files ────────────────────────────────────────────────────
    sol_file  = input("\nSolution file  (.h5): ").strip()
    mesh_file = input("Mesh file      (.h5): ").strip()
    out_file  = input("Output file    (.h5) [sol_modified.h5]: ").strip()
    if not out_file:
        out_file = "sol_modified.h5"

    if not os.path.isfile(sol_file):
        sys.exit(f"✗  Solution file not found: {sol_file}")
    if not os.path.isfile(mesh_file):
        sys.exit(f"✗  Mesh file not found: {mesh_file}")

    print(f"\nCopying {sol_file} → {out_file} …")
    shutil.copy2(sol_file, out_file)

    # ── Coordinates ──────────────────────────────────────────────
    print("\nLoading coordinates …")
    x, y, z = load_mesh_coords(mesh_file)
    print(f"  X: [{x.min():.4f} – {x.max():.4f}]")
    print(f"  Y: [{y.min():.4f} – {y.max():.4f}]")
    print(f"  Z: [{z.min():.4f} – {z.max():.4f}]")

    # ── Patch selection ───────────────────────────────────────────
    patch = select_patch()

    center  = patch["center"]
    axis    = patch["axis"]
    height  = patch["height"]
    r_inner = patch["r_inner"]
    r_outer = patch["r_outer"]
    thick   = patch["transition_thickness"]

    # ── Node selection ────────────────────────────────────────────
    print(f"\nSearching nodes in '{patch['name']}' …")
    idx_core = select_nodes_cylinder(x, y, z, center, axis, r_inner, r_outer, height)
    print(f"  → {len(idx_core):,} core nodes found")
    if len(idx_core) == 0:
        sys.exit("✗  No nodes found. Check coordinates and units.")

    print("Building transition zone …")
    idx_trans, alpha_trans = build_transition_zone(
        x, y, z, center, axis, r_inner, r_outer, height,
        transition_thickness=thick,
    )
    print(f"  → {len(idx_trans):,} transition nodes")

    # ── Species selection ─────────────────────────────────────────
    with h5py.File(out_file, "r") as f:
        species_list = list(f["RhoSpecies"].keys())
    total =0.0
    while total != 1.0:
        print("\n"+"──"*24)
        print("\n── Available species ──")
        for i, sp in enumerate(species_list, 1):
            print(f"  {i}. {sp}")

        species_changes = {}
        print("\nSpecies to modify (empty input = done):")
        while True:
            sp_input = input("  Species name: ").strip()
            if not sp_input:
                break
            if sp_input not in species_list:
                print(f"  ⚠  Valid species: {species_list}")
                continue
            val = prompt_float(f"  Y_{sp_input} = ")
            if not (0.0 <= val <= 1.0):
                print("  ⚠  Value must be between 0 and 1.")
                continue
            species_changes[sp_input] = val

        if not species_changes:
            sys.exit("✗  No species selected. Aborting.")

        total = sum(species_changes.values())
        free = [sp for sp in species_list if sp not in species_changes]
        print(f"\n  Imposed Y   : {species_changes}  (sum={total:.4f})")
        if total > 1.0:
            print(f"WARNING:  Sum of Y = {total:.4f} > 1. Please try again and make sure that SUM of Y = 1")
            total = 0.0
        if total < 1.0:
            print(f"WARNING:  Sum of Y = {total:.4f} < 1. Please try again and make sure that SUM of Y = 1")
    
    print(f"  Distributed over: {free}  (remaining={1-total:.4f})")

    # ── Summary ───────────────────────────────────────────────────
    cyl_str = "hollow" if r_inner > 0 else "solid"
    print("\n" + "─" * 55)
    print("  SUMMARY")
    print(f"  Patch          : {patch['name']}")
    print(f"  Output file    : {out_file}")
    print(f"  Center         : {center}")
    print(f"  Axis / Height  : {axis} / {height} m")
    print(f"  Cylinder       : {cyl_str}  r_inner={r_inner}  r_outer={r_outer}")
    print(f"  Transition     : {thick} m")
    print(f"  Core nodes     : {len(idx_core):,}")
    print(f"  Trans. nodes   : {len(idx_trans):,}")
    print(f"  Modifications  : {species_changes}")
    print("─" * 55)

    confirm = input("  Confirm modification? (yes/no): ").strip().lower()
    if confirm not in ("yes", "y"):
        sys.exit("Aborted.")

    # ── Apply ─────────────────────────────────────────────────────
    print("\nApplying modifications …")
    modify_species(out_file, idx_core, species_changes, idx_trans, alpha_trans)

    verify_solution(out_file, idx_core)

    print(f"\n HDF5 file saved: {out_file}")

    # ── XMF export ────────────────────────────────────────────────
    print("\n── XMF export for Paraview ──")
    xmf_mode = prompt_choice(
        "  Mode (dual = mesh+sol separate / single = all in sol): ",
        ["dual", "single"]
    )
    with h5py.File(out_file, "r") as f:
        t_phys = float(f["Parameters/dtsum"][0])
    xmf_out = os.path.splitext(out_file)[0] + ".xmf"
    generate_xmf(
        sol_file   = out_file,
        mesh_file  = mesh_file,
        xmf_file   = xmf_out,
        time_value = t_phys,
        mode       = xmf_mode,
    )
    print(f"\n XMF file saved: '{xmf_out}'")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    main()