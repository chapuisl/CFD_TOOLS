
"""
# ===================================================================================================================
#  CFD Solution Modifier - Patches cylinder for HYLON Geometry
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
#  This module stores all predefined patches predefined  used by the
#  CFD solution modifier.
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
# ═══════════════════════════════════════════════════════════════════
#  PREDEFINED PATCHES
#  Each patch is a dict with keys:
#    name, description,
#    center (x0,y0,z0), axis, height, r_inner, r_outer,
#    transition_thickness
# ═══════════════════════════════════════════════════════════════════

PREDEFINED_PATCHES = [
    {
        "name":        "HYLON_jet_central",
        "description": "HYLON central jet — solid cylinder",
        "center":      (0.0, 0.0, -0.1117),
        "axis":        "z",
        "height":      0.1067,
        "r_inner":     0.0,
        "r_outer":     0.003,
        "transition_thickness": 0.001,
    },
    {
        "name":        "HYLON_SWIRLER_ext",
        "description": "HYLON outer swirler — hollow cylinder",
        "center":      (0.0, 0.0, -0.1117),
        "axis":        "z",
        "height":      0.105,
        "r_inner":     0.0031,
        "r_outer":     0.04,
        "transition_thickness": 0.002,
    },
    {
        "name":        "combustion_chamber",
        "description": "Combustion chamber — solid cylinder",
        "center":      (0.0, 0.0, -0.002),
        "axis":        "z",
        "height":      0.15,
        "r_inner":     0.0,
        "r_outer":     0.05,
        "transition_thickness": 0.001,
    },
]
