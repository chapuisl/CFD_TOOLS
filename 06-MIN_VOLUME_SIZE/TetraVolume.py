"""
# ===================================================================================================================
#  Tetrahedron Volume Calculator - Main Program
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mecanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 03 September 2026
#  Last Modified   : 03 September 2026
#  Version         : 1.0.00
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  Compute the volume of a perfect (regular) tetrahedron from a given edge
#  length, and print a nicely formatted summary in the terminal showing:
#      - the reference volume V = sqrt(2)/12 * a^3
#      - the same volume expressed in m^3 and mm^3
#      - the volume obtained if the edge length is reduced by 5, 10, 15,
#        20 and 30 percent
#
# -------------------------------------------------------------------------------------------------------------------
#  USAGE EXAMPLES
# -------------------------------------------------------------------------------------------------------------------
#  python main.py 50            # edge length = 50 micrometers (default unit)
#  python main.py 50 um         # edge length = 50 micrometers (explicit)
#  python main.py 0.5 mm        # edge length = 0.5 millimeters
#  python main.py 0.0005 m        # edge length = 0.0005 meters
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
import time
import sys
import math
 
 
"""
# ===================================================================================================================
#  Timing decorator
# ===================================================================================================================
"""
 
timings = {}
 
#Description: Decorates a function to measure and record its execution time in the global timings dictionary.
#
#@func:   Function to time [callable]
#
#@return1: Wrapped function that times each call [callable]
def timer(func):
    def wrapper(*args, **kwargs):
        t0 = time.perf_counter()
        result = func(*args, **kwargs)
        timings[func.__name__] = time.perf_counter() - t0
        return result
    return wrapper
 
 
"""
# ===================================================================================================================
#  Configuration / constants
# ===================================================================================================================
"""
 
# Percentage reductions applied to the edge length (in %)
REDUCTIONS_PCT = [5, 10, 15, 20, 30]
 
# Conversion factor: 1 mm = 1000 um
MM_TO_UM = 1e3

# Conversion factor: 1 mm = 1000 um
M_TO_UM = 1e6
 
# Geometric constant for a regular tetrahedron: V = TETRA_CONST * a^3
TETRA_CONST = math.sqrt(2) / 12.0
 
 
"""
# ===================================================================================================================
#  Argument parsing
# ===================================================================================================================
"""
 
#Description: Parses command-line arguments (edge length and optional unit).
#
#@argv:   List of command-line arguments, excluding the program name [list[str]]
#
#@return1: Parsed arguments namespace with 'taille' [float] and 'unite' [str] [argparse.Namespace]
def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="Compute the volume of a perfect tetrahedron from an edge length."
    )
    parser.add_argument(
        "taille",
        type=float,
        help="Edge length of the tetrahedron (float)."
    )
    parser.add_argument(
        "unite",
        nargs="?",
        default="um",
        choices=["um", "mm","m"],
        help="Unit of the given edge length: 'um' (default, micrometers) or 'mm' (millimeters), or 'm' (meters)."
    )
    return parser.parse_args(argv)
 
 
"""
# ===================================================================================================================
#  Volume computation
# ===================================================================================================================
"""
 
#Description: Converts an edge length to micrometers, based on the given unit.
#
#@taille:   Edge length value [float]
#@unite:    Unit of the value, 'um' or 'mm' [str]
#
#@return1: Edge length expressed in micrometers [float]
def to_micrometers(taille, unite):
    if unite == "m":
        return taille * M_TO_UM

    if unite == "mm":
        return taille * MM_TO_UM

    else:
        return taille 

 
#Description: Computes the volume of a perfect tetrahedron of edge length a.
#
#@a:   Edge length [float]
#
#@return1: Volume of the tetrahedron, in the same unit as a, cubed [float]
@timer
def tetrahedron_volume(a):
    return TETRA_CONST * a ** 3
 
 
#Description: Builds the list of (percentage, reduced_edge, reduced_volume) for each
#             requested percentage reduction applied to the edge length.
#
#@a_um:          Reference edge length, in micrometers [float]
#@percentages:   List of percentage reductions to apply [list[int]]
#
#@return1: List of tuples (pct, reduced_edge_um, reduced_volume_um3) [list[tuple]]
def compute_reductions(a_um, percentages):
    results = []
    for pct in percentages:
        a_reduced = a_um * (1.0 - pct / 100.0)
        v_reduced = tetrahedron_volume(a_reduced)
        results.append((pct, a_reduced, v_reduced))
    return results
 
 
"""
# ===================================================================================================================
#  Pretty terminal output
# ===================================================================================================================
"""
 
#Description: Formats a volume expressed in um^3 into a human-readable string,
#             showing both um^3 and mm^3.
#
#@v_um3:   Volume in cubic micrometers [float]
#
#@return1: Formatted string "X.XXXXXXe+XX um^3  (Y.YYYYYYe+XX mm^3)" [str]
def format_volume(v_um3):
    v_mm3 = v_um3 / (MM_TO_UM ** 3)
    v_m3 = v_um3 / (M_TO_UM ** 3)
    return f"{v_m3:14.6e} m^3   |   {v_mm3:14.6e} mm^3"
 
 
#Description: Prints a nicely formatted summary of the tetrahedron volume and its
#             reductions to the terminal, using box-drawing characters.
#
#@a_input:    Original edge length, as given by the user [float]
#@unite:      Unit given by the user, 'um' or 'mm' [str]
#@a_um:       Edge length converted to micrometers [float]
#@v_um3:      Reference volume, in cubic micrometers [float]
#@reductions: List of (pct, reduced_edge_um, reduced_volume_um3) tuples [list[tuple]]
#
#@return1: None, console output only [None]
def print_summary(a_input, unite, a_um, v_um3, reductions):
    width = 78
    bar = "=" * width
 
    print(bar)
    print(" TETRAHEDRON VOLUME SUMMARY".center(width))
    print(bar)
    print(f" Edge length (input)   : {a_input:g} {unite}")
    print(f" Edge length (um)      : {a_um:.6f} um")
    print("-" * width)
    print(" Reference volume  V = sqrt(2)/12 * a^3")
    print(f"   {format_volume(v_um3)}")
    print(bar)
    print(" VOLUME FOR REDUCED EDGE LENGTHS".center(width))
    print(bar)
    header = f" {'Reduction':>10} | {'Edge length (um)':>18} | {'Volume':>34}"
    print(header)
    print("-" * width)
    for pct, a_reduced, v_reduced in reductions:
        print(f" {-pct:>9}% | {a_reduced:18.6f} | {format_volume(v_reduced)}")
    print(bar)
 
 
"""
# ===================================================================================================================
#  Main entry point
# ===================================================================================================================
"""
 
#Description: Program entry point: parses arguments, computes the tetrahedron volume
#             and its reductions, then prints the formatted summary.
#
#@argv:   List of command-line arguments, excluding the program name [list[str]]
#
#@return1: None [None]
def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
 
    if args.taille <= 0:
        print("Error: edge length must be strictly positive.", file=sys.stderr)
        sys.exit(1)
 
    a_um = to_micrometers(args.taille, args.unite)
    v_um3 = tetrahedron_volume(a_um)
    reductions = compute_reductions(a_um, REDUCTIONS_PCT)
 
    print_summary(args.taille, args.unite, a_um, v_um3, reductions)
 
 
if __name__ == "__main__":
    main()
