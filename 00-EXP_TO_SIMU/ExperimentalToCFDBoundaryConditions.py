"""
# ===================================================================================================================
#  Cantera Simulation Mechanism Loader
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mecanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 17 February 2026
#  Last Modified   : 03 July 2026
#  Version         : 2.0.00
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  Computes mass flow rates and mass fractions for an HYLON injection configuration
#  (central H2 jet + NH3/air coflow), from two possible input types:
#
#    - Type "power" : P, T, PhiGlobal, alpha_H2 (share of thermal power supplied by H2), Thermic_Power
#    - Type "NL"     : P, T, m_AIR_NL, m_NH3_NL, m_H2_NL (flow rates in Nl/min)
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
#  This configuration file is intended solely for use within the associated
#  simulation framework developed by the author.
#
# ===================================================================================================================
"""

import cantera as ct
import numpy as np

"""
# ===================================================================================================================
#  Base functions (unchanged logic)
# ===================================================================================================================
"""

#DESCRIPTION: Computes the volume occupied by 1 mole of a perfect gas at a given pressure and temperature.
#
#@input  P: pressure [Pa]
#@input  T: temperature [K]
#
#@return Vperfect: molar volume [L]
def PerfectVolume_PT_for_1mol(P, T):
    Vperfect = 8.314 * T / P * 1000  # L
    return Vperfect


#DESCRIPTION: Converts a volumetric flow rate given in normal liters per minute (Nl/min) into a mass
#             flow rate in kg/s, using the standard density of the species at 0 degC / 1 atm (source: Air Liquide).
#
#@input  Flow_rate_NL: flow rate [Nl/min]
#@input  species: species name, one of 'air', 'H2', 'CH4', 'NH3'
#
#@return debit_kgs: mass flow rate [kg/s]
def NormoLitre_TO_kg(Flow_rate_NL, species):
    # densities at 0 degC and 1 atm - source: Air Liquide
    rho0_air = 1.292   # g/L
    rho0_CH4 = 0.7173  # g/L
    rho0_NH3 = 0.7713  # g/L
    rho0_H2 = 0.0899   # g/L
    gas_properties = {
        'air': {'rho': rho0_air},
        'H2':  {'rho': rho0_H2},
        'CH4': {'rho': rho0_CH4},
        'NH3': {'rho': rho0_NH3},
    }
    rho = gas_properties[species]['rho']
    debit_kgs = (Flow_rate_NL * rho) / (60 * 1000)
    return debit_kgs


#DESCRIPTION: Converts a mass flow rate given in kg/s into a volumetric flow rate in normal liters per
#             minute (Nl/min), using the standard density of the species at 0 degC / 1 atm (source: Air Liquide).
#
#@input  Flow_rate_Kg: mass flow rate [kg/s]
#@input  species: species name, one of 'air', 'H2', 'CH4', 'NH3'
#
#@return debit_NLmin: flow rate [Nl/min]
def Kg_TO_NormoLitre(Flow_rate_Kg, species):
    # densities at 0 degC and 1 atm - source: Air Liquide
    rho0_air = 1.292   # g/L
    rho0_CH4 = 0.7173  # g/L
    rho0_NH3 = 0.7713  # g/L
    rho0_H2 = 0.0899   # g/L
    gas_properties = {
        'air': {'rho': rho0_air},
        'H2':  {'rho': rho0_H2},
        'CH4': {'rho': rho0_CH4},
        'NH3': {'rho': rho0_NH3},
    }
    rho = gas_properties[species]['rho']
    debit_NLmin = (Flow_rate_Kg * 60 * 1000) / rho
    return debit_NLmin


#DESCRIPTION: Computes the oxygen mass fraction in air, assuming a simple O2/N2 composition.
#
#@return Y_O2: oxygen mass fraction in air [-]
def Y_O2air():
    W_O2 = 0.032  # kg/mol
    W_N2 = 0.028  # kg/mol
    n_O2 = 0.21   # mol
    n_N2 = 1 - n_O2  # mol
    Y_O2 = (n_O2 * W_O2) / (n_O2 * W_O2 + n_N2 * W_N2)
    return Y_O2


#DESCRIPTION: Computes the stoichiometric air-to-fuel mass ratio of a given fuel (or fuel blend) against
#             pure O2, using the GRI-3.0 mechanism in Cantera.
#
#@input  fuel: Cantera fuel composition string (e.g. "NH3:1.0" or "NH3:0.3,H2:0.7")
#
#@return s: stoichiometric air-to-fuel mass ratio [-]
def stoich_rapport(fuel):
    gas = ct.Solution("gri30.yaml")
    s = gas.stoich_air_fuel_ratio(fuel=fuel, oxidizer="O2:1")
    return s


#DESCRIPTION: Computes the fuel mass flow rate required to reach a target equivalence ratio Phi in a
#             premixed air/fuel stream, given the air mass flow rate.
#
#@input  m_AIR: air mass flow rate [kg/s]
#@input  Phi: target equivalence ratio [-]
#@input  fuel: Cantera fuel composition string (e.g. "NH3:1.0")
#
#@return m_F: fuel mass flow rate [kg/s]
def massAIR_To_massFuel_in_premixed(m_AIR, Phi, fuel):
    s = stoich_rapport(fuel)
    Y_O2 = Y_O2air()
    m_F = Phi / s * m_AIR * Y_O2
    return m_F


#DESCRIPTION: Computes the fuel mass flow rate corresponding to a given thermal power output, using the
#             fuel's lower heating value (LHV / PCI).
#
#@input  Power: thermal power [kW]
#@input  PCI: lower heating value of the fuel [MJ/kg]
#
#@return debit: fuel mass flow rate [kg/s]
def massF_FROM_Thermic_POWER(Power, PCI):
    Power = Power * 1000
    PCI = PCI * 1e6
    return Power / PCI


"""
# ===================================================================================================================
#  Conversion of both input types into a common set of mass flow rates (kg/s)
# ===================================================================================================================
"""

#DESCRIPTION: Input type "power". Computes air, H2 and NH3 mass flow rates from the total thermal power,
#             the share of power supplied by H2 (alpha_H2) and the targeted global equivalence ratio.
#             Method: (1) split the thermal power between H2 and NH3 using alpha_H2 ; (2) get m_H2 and
#             m_NH3 from their respective LHV ; (3) build the NH3/H2 blend and get its stoichiometric
#             ratio s ; (4) deduce m_AIR from PhiGlobal = s*(m_H2+m_NH3)/(m_AIR*Y_O2).
#
#@input  Thermic_Power: total thermal power [kW]
#@input  alpha_H2: share of thermal power supplied by H2 [-] (between 0 and 1)
#@input  PhiGlobal: targeted global equivalence ratio [-]
#@input  PCI_H2: lower heating value of H2 [MJ/kg]
#@input  PCI_NH3: lower heating value of NH3 [MJ/kg]
#
#@return m_AIR: air mass flow rate [kg/s]
#@return m_H2: H2 mass flow rate [kg/s]
#@return m_NH3: NH3 mass flow rate [kg/s]
def flows_from_thermal_power(Thermic_Power, alpha_H2, PhiGlobal, PCI_H2, PCI_NH3):
    P_H2 = alpha_H2 * Thermic_Power
    P_NH3 = Thermic_Power - P_H2

    m_H2 = massF_FROM_Thermic_POWER(P_H2, PCI_H2)
    m_NH3 = massF_FROM_Thermic_POWER(P_NH3, PCI_NH3)

    W_H2, W_NH3 = 0.002, 0.017  # kg/mol
    Beta = (m_H2 / W_H2) / ((m_H2 / W_H2) + (m_NH3 / W_NH3))
    fuel_blend = f"NH3:{1 - Beta},H2:{Beta}"
    s_blend = stoich_rapport(fuel_blend)
    Y_O2 = Y_O2air()

    m_AIR = s_blend * (m_H2 + m_NH3) / (PhiGlobal * Y_O2)

    return m_AIR, m_H2, m_NH3


#DESCRIPTION: Input type "NL". Converts air, H2 and NH3 flow rates given in Nl/min into mass flow rates in kg/s.
#
#@input  m_AIR_NL: air flow rate [Nl/min]
#@input  m_H2_NL: H2 flow rate [Nl/min]
#@input  m_NH3_NL: NH3 flow rate [Nl/min]
#
#@return m_AIR: air mass flow rate [kg/s]
#@return m_H2: H2 mass flow rate [kg/s]
#@return m_NH3: NH3 mass flow rate [kg/s]
def flows_from_NL(m_AIR_NL, m_H2_NL, m_NH3_NL):
    m_AIR = NormoLitre_TO_kg(m_AIR_NL, 'air')
    m_H2 = NormoLitre_TO_kg(m_H2_NL, 'H2')
    m_NH3 = NormoLitre_TO_kg(m_NH3_NL, 'NH3')
    return m_AIR, m_H2, m_NH3


"""
# ===================================================================================================================
#  Common computation of the HYLON injection state (mass fractions, equivalence ratios) from mass flow rates
# ===================================================================================================================
"""

#DESCRIPTION: From the air, H2 and NH3 mass flow rates, computes the fuel blend composition, the global
#             equivalence ratio, the mass fractions in the central jet (pure H2) and in the coflow
#             (NH3 + air), and the coflow equivalence ratio.
#
#@input  m_AIR: air mass flow rate [kg/s]
#@input  m_H2: H2 mass flow rate [kg/s]
#@input  m_NH3: NH3 mass flow rate [kg/s]
#
#@return state: dict containing Beta, fuel_blend, PhiGlobal, PhiCoflow, Y_jet, Y_coflow, m_AIR, m_H2, m_NH3
def compute_HYLON_state(m_AIR, m_H2, m_NH3):
    W_H2, W_NH3 = 0.002, 0.017  # kg/mol
    Beta = (m_H2 / W_H2) / ((m_H2 / W_H2) + (m_NH3 / W_NH3))
    fuel_blend = f"NH3:{1 - Beta},H2:{Beta}"
    s_blend = stoich_rapport(fuel_blend)
    Y_O2 = Y_O2air()

    PhiGlobal_calc = s_blend * (m_H2 + m_NH3) / (m_AIR * Y_O2)

    # Central jet: pure H2
    Y_jet = {"H2": 1.0, "NH3": 0.0, "N2": 0.0, "O2": 0.0}

    # Coflow: NH3 + air
    Y_NH3_cf = m_NH3 / (m_NH3 + m_AIR)
    Y_O2_cf = (m_AIR * Y_O2) / (m_NH3 + m_AIR)
    Y_N2_cf = (m_AIR * (1 - Y_O2)) / (m_NH3 + m_AIR)
    Y_coflow = {"H2": 0.0, "NH3": Y_NH3_cf, "N2": Y_N2_cf, "O2": Y_O2_cf}

    s_NH3 = stoich_rapport("NH3:1.0")
    PhiCoflow_calc = s_NH3 * Y_NH3_cf / Y_O2_cf

    return {
        "Beta": Beta,
        "fuel_blend": fuel_blend,
        "PhiGlobal": PhiGlobal_calc,
        "PhiCoflow": PhiCoflow_calc,
        "Y_jet": Y_jet,
        "Y_coflow": Y_coflow,
        "m_AIR": m_AIR,
        "m_H2": m_H2,
        "m_NH3": m_NH3,
    }


#DESCRIPTION: Prints a formatted summary of an HYLON injection state (mass flow rates, fuel blend,
#             global equivalence ratio, jet and coflow mass fractions, coflow equivalence ratio),
#             optionally comparing the computed equivalence ratios to experimental reference values.
#
#@input  state: dict returned by compute_HYLON_state
#@input  PhiGlobal_ref: experimental reference value for the global equivalence ratio [-] (optional, default None)
#@input  PhiCoflow_ref: experimental reference value for the coflow equivalence ratio [-] (optional, default None)
#
#@return None
def print_HYLON_state(state, PhiGlobal_ref=None, PhiCoflow_ref=None):
    m_AIR, m_H2, m_NH3 = state["m_AIR"], state["m_H2"], state["m_NH3"]

    print("\n" + "=" * 70)
    print("                         MASS FLOW RATES")
    print("=" * 70)
    print(f"{'Species':<8} | {'mdot [g/s]':>12} | {'mdot [kg/s]':>12} | {'mdot [NL/min]':>15}")
    print("-" * 70)
    print(f"{'H2':<8} | {m_H2 * 1000:12.3f} | {m_H2 :12.3e}| {Kg_TO_NormoLitre(m_H2, 'H2'):15.2f}")
    print(f"{'NH3':<8} | {m_NH3 * 1000:12.3f} | {m_NH3 :12.3e}| {Kg_TO_NormoLitre(m_NH3, 'NH3'):15.2f}")
    print(f"{'Air':<8} | {m_AIR * 1000:12.3f} | {m_AIR :12.3e}| {Kg_TO_NormoLitre(m_AIR, 'air'):15.2f}")
    print(f"{'Air+ NH3':<8} | {(m_AIR +m_NH3) * 1000:12.3f} | {m_AIR + m_NH3 :12.3e}| {Kg_TO_NormoLitre(m_AIR, 'air')+Kg_TO_NormoLitre(m_NH3, 'NH3'):15.2f}")
    print("=" * 70)

    print(f"\nFuel blend composition: {state['fuel_blend']}")
    ref_txt = f"  (exp reference: {PhiGlobal_ref})" if PhiGlobal_ref is not None else ""
    print(f"Computed global equivalence ratio: {state['PhiGlobal']:.4f}{ref_txt}")

    print("\n--- MASS FRACTIONS - Central jet (H2) ---")
    for k, v in state["Y_jet"].items():
        print(f"  * {k:<4} : {v:10.4f}")

    print("\n--- MASS FRACTIONS - Coflow (NH3 + air) ---")
    for k, v in state["Y_coflow"].items():
        print(f"  * {k:<4} : {v:10.4f}")

    ref_txt = f"  (exp reference: {PhiCoflow_ref})" if PhiCoflow_ref is not None else ""
    print(f"\nComputed coflow equivalence ratio: {state['PhiCoflow']:.4f}{ref_txt}")


"""
# ===================================================================================================================
#  Single entry point: accepts both input types
# ===================================================================================================================
"""

#DESCRIPTION: Single entry point for the HYLON injection calculation. Dispatches to the appropriate
#             conversion function depending on input_type, computes the injection state (mass flow
#             rates, mass fractions, equivalence ratios) and prints a formatted summary.
#
#@input  input_type: "power" or "NL"
#@input  kw: keyword arguments, depending on input_type -
#            "power" requires Thermic_Power [kW], alpha_H2 [-], PhiGlobal [-], PCI_H2 [MJ/kg], PCI_NH3 [MJ/kg],
#                     and optionally accepts PhiCoflow [-] as an experimental reference value ;
#            "NL" requires m_AIR_NL, m_H2_NL, m_NH3_NL [Nl/min],
#                 and optionally accepts PhiGlobal [-] and PhiCoflow [-] as experimental reference values.
#
#@return state: dict returned by compute_HYLON_state
def Run_HYLON_Injection(input_type, **kw):
    if input_type == "power":
        m_AIR, m_H2, m_NH3 = flows_from_thermal_power(
            kw["Thermic_Power"], kw["alpha_H2"], kw["PhiGlobal"],
            kw["PCI_H2"], kw["PCI_NH3"]
        )
        PhiGlobal_ref = kw["PhiGlobal"]
        PhiCoflow_ref = kw.get("PhiCoflow")

    elif input_type == "NL":
        m_AIR, m_H2, m_NH3 = flows_from_NL(
            kw["m_AIR_NL"], kw["m_H2_NL"], kw["m_NH3_NL"]
        )
        PhiGlobal_ref = kw.get("PhiGlobal")
        PhiCoflow_ref = kw.get("PhiCoflow")

    else:
        raise ValueError('input_type must be "power" or "NL"')

    state = compute_HYLON_state(m_AIR, m_H2, m_NH3)
    print_HYLON_state(state, PhiGlobal_ref, PhiCoflow_ref)
    return state


"""
# ===================================================================================================================
#  Usage example
# ===================================================================================================================
"""

def __main__():
     print("########## CASE 1: Power / Equivalence ratio inputs (P, T, PhiGlobal, alpha_H2, Pth) ##########")
     state1 = Run_HYLON_Injection(
        "power",
        Thermic_Power=6,     # kW
        alpha_H2=0.7,
        PhiGlobal=0.5,
        PCI_H2=120,          # MJ/kg
        PCI_NH3=18.6,        # MJ/kg
        PhiCoflow=0.096,     # experimental reference value (optional, for comparison)
    )
#      print("\n\n########## CASE 2: Flow rate inputs (P, T, m_AIR_NL, m_NH3_NL, m_H2_NL) ##########")
#      state2 = Run_HYLON_Injection(
#         "NL",
#         m_AIR_NL=378,
#         m_H2_NL=31.15,
#         m_NH3_NL=10.04,
#     )


__main__()