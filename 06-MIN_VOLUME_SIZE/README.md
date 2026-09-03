# Tetrahedron Volume Calculator

A small command-line tool that computes the volume of a **perfect (regular) tetrahedron**
from a given edge length, and prints a nicely formatted summary in the terminal —
including the volume you'd get if the edge length were reduced by 5, 10, 15, 20 or 30 %.

## Purpose

Given an edge length `a`, the script computes the reference volume of a regular
tetrahedron using the exact geometric formula:

```
V = (sqrt(2) / 12) * a^3
```

It then reports:

- the edge length you entered, converted to micrometers (µm),
- the reference volume, expressed in both µm³ and mm³,
- the volume obtained for the same tetrahedron shape but with the edge length
  reduced by 5 %, 10 %, 15 %, 20 % and 30 % (useful, for example, to see how much
  smaller a mesh cell's tetrahedron becomes as its characteristic size shrinks).

## Requirements

- Python 3.6 or later
- No external dependencies — only the Python standard library (`argparse`, `math`,
  `time`, `sys`) is used.

## Usage

```bash
python main.py <taille> [unite]
```

| Argument | Required | Description                                              |
|----------|----------|-----------------------------------------------------------|
| `taille` | Yes      | Edge length of the tetrahedron (a positive float).         |
| `unite`  | No       | Unit of `taille`: `um` (micrometers, **default**) or `mm`. |

### Examples

Edge length of 50 micrometers (default unit, no need to specify `um`):

```bash
python main.py 50
```

Edge length of 50 micrometers (unit given explicitly):

```bash
python main.py 50 um
```

Edge length of 0.5 millimeters:

```bash
python main.py 0.5 mm
```

### Sample output

```
==============================================================================
                          TETRAHEDRON VOLUME SUMMARY
==============================================================================
 Edge length (input)   : 50 um
 Edge length (um)      : 50.000000 um
------------------------------------------------------------------------------
 Reference volume  V = sqrt(2)/12 * a^3
     1.473139e+04 um^3   |     1.473139e-05 mm^3
==============================================================================
                        VOLUME FOR REDUCED EDGE LENGTHS
==============================================================================
  Reduction |   Edge length (um) |                             Volume
------------------------------------------------------------------------------
        -5% |          47.500000 |   1.263033e+04 um^3   |     1.263033e-05 mm^3
       -10% |          45.000000 |   1.073918e+04 um^3   |     1.073918e-05 mm^3
       -15% |          42.500000 |   9.046916e+03 um^3   |     9.046916e-06 mm^3
       -20% |          40.000000 |   7.542472e+03 um^3   |     7.542472e-06 mm^3
       -30% |          35.000000 |   5.052867e+03 um^3   |     5.052867e-06 mm^3
==============================================================================
```

## How it works

The script is organized into a few clear stages:

1. **Argument parsing** (`parse_args`) — reads the edge length and the optional
   unit (`um` or `mm`) from the command line using `argparse`.
2. **Unit conversion** (`to_micrometers`) — converts the input edge length to
   micrometers, the internal working unit.
3. **Volume computation** (`tetrahedron_volume`) — applies the exact formula
   `V = sqrt(2)/12 * a^3` to compute the volume of a regular tetrahedron.
4. **Reductions** (`compute_reductions`) — for each percentage in `[5, 10, 15, 20, 30]`,
   shrinks the edge length by that percentage and recomputes the resulting volume.
5. **Formatted output** (`format_volume`, `print_summary`) — prints a clean,
   boxed summary table to the terminal, showing every volume in both µm³ and mm³.

Each core computation function is wrapped with a `@timer` decorator that records
its execution time in the global `timings` dictionary, in case you want to inspect
performance later (e.g., by printing `timings` after calling `main()`).

## Notes

- The edge-length reductions are applied to the **edge length**, not directly to
  the volume — since volume scales with the cube of the edge length, a 10 %
  reduction in edge length results in roughly a 27 % reduction in volume, not 10 %.
- The script exits with an error message if a non-positive edge length is provided.
