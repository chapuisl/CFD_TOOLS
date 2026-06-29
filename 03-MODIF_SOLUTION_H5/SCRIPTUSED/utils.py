"""
# ===================================================================================================================
#  CFD Solution Modifier - Utils function
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
#  This module stores all little function usefull for the CFD solution modifier.
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
import h5py

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def prompt_float(msg):
    while True:
        raw = input(msg).strip()
        try:
            return float(raw)
        except ValueError:
            print("  ⚠  A number is expected.")


def prompt_float_opt(msg, default):
    """Prompt with default value (Enter = keep current)."""
    while True:
        raw = input(msg).strip()
        if raw == "":
            return default
        try:
            return float(raw)
        except ValueError:
            print("  ⚠  A number is expected.")


def prompt_choice(msg, choices):
    while True:
        raw = input(msg).strip().lower()
        if raw in choices:
            return raw
        print(f"  ⚠  Valid choices: {choices}")


def verify_solution(output_file, idx_core):
    """Print a verification summary after modification."""
    print("\n─── Verification ───")
    with h5py.File(output_file, "r") as f:
        rho     = f["GaseousPhase/rho"][:]
        species = list(f["RhoSpecies"].keys())
        sum_Y   = sum(f[f"RhoSpecies/{sp}"][:] for sp in species) / rho
        print(f"  sum(Y_k) global : [{sum_Y.min():.6f} – {sum_Y.max():.6f}]  (target ~1)")
        print(f"  sum(Y_k) core   : [{sum_Y[idx_core].min():.6f} – {sum_Y[idx_core].max():.6f}]")
        print(f"  rho global      : [{rho.min():.4f} – {rho.max():.4f}] kg/m³")
        print(f"  rho core        : [{rho[idx_core].min():.4f} – {rho[idx_core].max():.4f}] kg/m³")
        print("\n  Core mass fractions:")
        for sp in species:
            Y_c = f[f"RhoSpecies/{sp}"][:][idx_core] / rho[idx_core]
            print(f"    Y_{sp:<6s}: [{Y_c.min():.6f} – {Y_c.max():.6f}]  mean={Y_c.mean():.6f}")

