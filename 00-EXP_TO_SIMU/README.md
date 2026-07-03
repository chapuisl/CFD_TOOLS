# HYLON Injection

Computation of mass flow rates and mass fractions for the **HYLON** injection configuration
(central H2 jet + NH3/air coflow), based on [Cantera](https://cantera.org/).

The script accepts two different sets of input parameters and returns the same output: mass
flow rates (kg/s and Nl/min), fuel blend composition, jet/coflow mass fractions, and global /
coflow equivalence ratios.

## Table of contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration overview](#configuration-overview)
- [Code architecture](#code-architecture)
- [Usage](#usage)
- [Function reference](#function-reference)
- [Author](#author)

## Requirements

| Library   | Tested version | Purpose                                   |
|-----------|----------------|--------------------------------------------|
| Python    | ≥ 3.9          | Language runtime                          |
| [Cantera](https://cantera.org/) | ≥ 3.2 | Chemical kinetics (stoichiometric ratios, `gri30.yaml` mechanism) |
| [NumPy](https://numpy.org/)     | ≥ 2.0 | Numerical utilities |

## Installation

```bash
pip install cantera numpy
```

> Cantera ships the `gri30.yaml` mechanism used by `stoich_rapport()`, so no extra
> mechanism file needs to be downloaded.

## Configuration overview

HYLON is a burner configuration with two concentric streams:

- **Central jet** — pure H2.
- **Coflow** — a premixed NH3/air stream.

The two streams share the same fuel blend on a *global* basis (jet + coflow combined), which
defines the **global equivalence ratio** `PhiGlobal`. The coflow alone, being leaner (diluted
by air), has its own **coflow equivalence ratio** `PhiCoflow`.

## Code architecture

The script is organized in four layers, so that both input types converge to a single
calculation path:

```
┌─────────────────────────────┐   ┌─────────────────────────────┐
│   Input type "power"        │   │   Input type "NL"           │
│   (Pth, alpha_H2, PhiGlobal)│   │   (m_AIR_NL, m_H2_NL,        │
│                              │   │    m_NH3_NL)                 │
└──────────────┬───────────────┘   └──────────────┬───────────────┘
               │ flows_from_thermal_power()        │ flows_from_NL()
               └──────────────────┬─────────────────┘
                                   ▼
                     common mass flow rates (kg/s)
                        m_AIR, m_H2, m_NH3
                                   │
                                   ▼
                       compute_HYLON_state()
              (fuel blend, PhiGlobal, PhiCoflow,
                     Y_jet, Y_coflow)
                                   │
                                   ▼
                       print_HYLON_state()
                    (formatted console output)
```

1. **Base functions** — unit conversions and thermochemistry helpers
   (`NormoLitre_TO_kg`, `Kg_TO_NormoLitre`, `Y_O2air`, `stoich_rapport`,
   `massAIR_To_massFuel_in_premixed`, `massF_FROM_Thermic_POWER`, `PerfectVolume_PT_for_1mol`).
2. **Input converters** — turn either input type into a common set of mass flow rates
   (`flows_from_thermal_power`, `flows_from_NL`).
3. **State calculation** — from the mass flow rates, computes the fuel blend, the mass
   fractions of each stream, and both equivalence ratios (`compute_HYLON_state`).
4. **Entry point / reporting** — dispatches to the right converter and prints a formatted
   summary (`Run_HYLON_Injection`, `print_HYLON_state`).

## Usage

### Case 1 — Power / equivalence ratio inputs

Given the total thermal power, the share of that power supplied by H2, and the targeted
global equivalence ratio:

```python
from HYLON_Injection import Run_HYLON_Injection

state = Run_HYLON_Injection(
    "power",
    Thermic_Power=8,     # kW
    alpha_H2=0.7,        # share of power from H2
    PhiGlobal=0.29,
    PCI_H2=120,          # MJ/kg
    PCI_NH3=18.6,         # MJ/kg
    PhiCoflow=0.096,      # optional experimental reference, printed for comparison
)
```

### Case 2 — Flow rate inputs

Given the air, H2 and NH3 flow rates directly, in Nl/min:

```python
state = Run_HYLON_Injection(
    "NL",
    m_AIR_NL=378,
    m_H2_NL=31.15,
    m_NH3_NL=10.04,
)
```

Both calls return a `state` dictionary:

```python
{
    "Beta": ...,          # H2 mole share within the fuel blend
    "fuel_blend": ...,    # Cantera fuel string, e.g. "NH3:0.25,H2:0.75"
    "PhiGlobal": ...,     # computed global equivalence ratio
    "PhiCoflow": ...,     # computed coflow equivalence ratio
    "Y_jet": {...},       # mass fractions in the central jet
    "Y_coflow": {...},    # mass fractions in the coflow
    "m_AIR": ...,         # air mass flow rate [kg/s]
    "m_H2": ...,          # H2 mass flow rate [kg/s]
    "m_NH3": ...,         # NH3 mass flow rate [kg/s]
}
```

Run the file directly to see both cases printed to the console:

```bash
python3 HYLON_Injection.py
```

## Function reference

| Function | Role |
|---|---|
| `PerfectVolume_PT_for_1mol(P, T)` | Molar volume of an ideal gas at `(P, T)` |
| `NormoLitre_TO_kg(Flow_rate_NL, species)` | Nl/min → kg/s |
| `Kg_TO_NormoLitre(Flow_rate_Kg, species)` | kg/s → Nl/min |
| `Y_O2air()` | O2 mass fraction in air |
| `stoich_rapport(fuel)` | Stoichiometric air/fuel mass ratio (Cantera, `gri30.yaml`) |
| `massAIR_To_massFuel_in_premixed(m_AIR, Phi, fuel)` | Fuel flow rate for a target `Phi` in a premixed stream |
| `massF_FROM_Thermic_POWER(Power, PCI)` | Fuel flow rate from thermal power and LHV |
| `flows_from_thermal_power(...)` | Input type "power" → `(m_AIR, m_H2, m_NH3)` |
| `flows_from_NL(...)` | Input type "NL" → `(m_AIR, m_H2, m_NH3)` |
| `compute_HYLON_state(m_AIR, m_H2, m_NH3)` | Mass flow rates → full injection state |
| `print_HYLON_state(state, ...)` | Formatted console report |
| `Run_HYLON_Injection(input_type, **kw)` | Single entry point (dispatch + compute + print) |

## Author

**Lilian CHAPUIS** — IMFT (Institut de Mécanique des Fluides de Toulouse), Toulouse, France.

© 2026 Lilian CHAPUIS — All Rights Reserved. See header of `HYLON_Injection.py` for the
copyright notice.
