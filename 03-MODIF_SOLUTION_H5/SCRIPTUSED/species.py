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
import numpy as np
import h5py
"""
# ===================================================================================================================
#  Imports from other modules
# ===================================================================================================================
"""

from SCRIPTUSED.constant import R_UNIVERSAL, MOLAR_MASSES


"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""


def compute_r_bar(Y, species_list):
    """
    r_bar = R_univ * sum(Y_k / M_k)   [J/kg/K]
    Mixture-specific gas constant.
    """
    inv_M_mix = np.zeros_like(list(Y.values())[0], dtype=np.float64)
    for sp in species_list:
        if sp not in MOLAR_MASSES:
            raise KeyError(
                f"Unknown molar mass for '{sp}'. "
                f"Add it to the MOLAR_MASSES dictionary."
            )
        inv_M_mix += Y[sp] / MOLAR_MASSES[sp]
    return R_UNIVERSAL * inv_M_mix


def modify_species(output_file, idx_core, species_changes, idx_trans, alpha_trans):
    """
    Apply new mass fractions to core nodes, apply a linear continuous
    gradient over the transition zone, then update all conservative variables.

    Parameters
    ----------
    output_file    : path to the HDF5 output file (writable copy)
    idx_core       : node indices inside the cylinder core
    species_changes: {species_name: new_Y_value}
                     unmodified species are renormalised automatically
    idx_trans      : transition zone node indices
    alpha_trans    : per-node linear weight in [0, 1]
    """
    with h5py.File(output_file, "r+") as f:
        species_list = list(f["RhoSpecies"].keys())

        # Read all current mass fractions
        rho_orig = f["GaseousPhase/rho"][:]
        Y_orig   = {sp: f[f"RhoSpecies/{sp}"][:] / rho_orig for sp in species_list}
        Y_new    = {sp: Y_orig[sp].copy() for sp in species_list}

        fixed_species = list(species_changes.keys())
        fixed_sum     = sum(species_changes.values())
        if fixed_sum > 1.0 + 1e-10:
            raise ValueError(f"Sum of imposed mass fractions = {fixed_sum:.4f} > 1!")

        free_species     = [sp for sp in species_list if sp not in fixed_species]
        free_sum_at_core = sum(Y_orig[sp][idx_core] for sp in free_species)
        remaining        = 1.0 - fixed_sum

        # ── Core: impose new values ───────────────────────────────
        for sp in fixed_species:
            Y_new[sp][idx_core] = species_changes[sp]

        # Renormalise free species proportionally
        for sp in free_species:
            denom = free_sum_at_core
            ratio = np.where(
                denom > 0,
                Y_orig[sp][idx_core] / denom * remaining,
                remaining / max(len(free_species), 1)
            )
            Y_new[sp][idx_core] = ratio

        # ── Transition: continuous linear interpolation per node ──
        if len(idx_trans) > 0:
            # Target value at alpha=0 (core boundary side)
            Y_inner = {sp: float(species_changes[sp]) for sp in fixed_species}
            free_sum_g = sum(float(Y_orig[sp][idx_core].mean()) for sp in free_species)
            for sp in free_species:
                mean_sp = float(Y_orig[sp][idx_core].mean())
                Y_inner[sp] = (mean_sp / free_sum_g * remaining
                               if free_sum_g > 0
                               else remaining / max(len(free_species), 1))

            for sp in species_list:
                # alpha=0 → Y_inner (core boundary)
                # alpha=1 → Y_orig  (unmodified exterior)
                Y_new[sp][idx_trans] = (
                    (1.0 - alpha_trans) * Y_inner[sp]
                    + alpha_trans       * Y_orig[sp][idx_trans]
                )

        # ── Combine all modified nodes ────────────────────────────
        all_idx = np.unique(
            np.concatenate([idx_core, idx_trans]) if len(idx_trans) > 0 else idx_core
        )

        # ── Update conservative variables ─────────────────────────
        rho  = f["GaseousPhase/rho"][:]
        rhoE = f["GaseousPhase/rhoE"][:]
        rhou = f["GaseousPhase/rhou"][:]
        rhov = f["GaseousPhase/rhov"][:]
        rhow = f["GaseousPhase/rhow"][:]
        P    = f["Additionals/pressure"][:].astype(np.float64)
        T    = f["Additionals/temperature"][:].astype(np.float64)

        # Store total specific energy and velocity BEFORE modification
        e_tot_old = rhoE[all_idx] / rho[all_idx]
        u_old     = rhou[all_idx] / rho[all_idx]
        v_old     = rhov[all_idx] / rho[all_idx]
        w_old     = rhow[all_idx] / rho[all_idx]

        # New rho via ideal gas law: rho = P / (r_bar * T)
        Y_idx     = {sp: Y_new[sp][all_idx] for sp in species_list}
        r_bar_new = compute_r_bar(Y_idx, species_list)
        rho_new   = P[all_idx] / (r_bar_new * T[all_idx])

        # Write new rhoY_k
        for sp in species_list:
            rhoY = f[f"RhoSpecies/{sp}"][:]
            rhoY[all_idx] = Y_new[sp][all_idx] * rho_new
            f[f"RhoSpecies/{sp}"][:] = rhoY

        # Write new rho, rhoE (e_tot conserved), rhou/rhov/rhow (velocity conserved)
        rho[all_idx]  = rho_new
        rhoE[all_idx] = rho_new * e_tot_old
        rhou[all_idx] = rho_new * u_old
        rhov[all_idx] = rho_new * v_old
        rhow[all_idx] = rho_new * w_old

        f["GaseousPhase/rho"][:] = rho
        f["GaseousPhase/rhoE"][:] = rhoE
        f["GaseousPhase/rhou"][:] = rhou
        f["GaseousPhase/rhov"][:] = rhov
        f["GaseousPhase/rhow"][:] = rhow

        # Update r_bar in Additionals if present
        if "r_bar" in f["Additionals"]:
            r_bar_arr = f["Additionals/r_bar"][:]
            r_bar_arr[all_idx] = r_bar_new.astype(np.float32)
            f["Additionals/r_bar"][:] = r_bar_arr

        print(f"\n  ✓ Core nodes modified      : {len(idx_core):,}")
        print(f"  ✓ Transition nodes modified : {len(idx_trans):,}")
        print(f"  ✓ Total nodes updated       : {len(all_idx):,}")
        print(f"  ✓ rho range                 : [{rho_new.min():.4f} – {rho_new.max():.4f}] kg/m³")
