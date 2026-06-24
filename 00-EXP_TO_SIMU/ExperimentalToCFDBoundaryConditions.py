"""
# ===================================================================================================================
#  Cantera Simulation Mechanism Loader
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mécanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 17 February 2026
#  Last Modified   : 18 February 2026
#  Version         : 1.0.01
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  This function provides a centralized registry of kinetic mechanisms that can be used in Cantera-based
#  combust
# -------------------------------------------------------------------------------------------------------------------
#  COPYRIGHT NOTICE
# -------------------------------------------------------------------------------------------------------------------
#  © 2026 Lilian CHAPUIS – All Rights Reserved.
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
"""
# ===================================================================================================================
#  Import library
# ===================================================================================================================
"""

import cantera as ct
import numpy as np

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def PerfectVolume_PT_for_1mol(P,T):
     Vperfect = 8.314 *T/P *1000 #L
     return Vperfect


def NormoLitre_TO_kg(Flow_rate_NL,species):
     # masse volumique à 0°C et 1 atm — source : Air Liquide
     rho0_air = 1.292  #g/L
     rho0_CH4 = 0.7173 #g/L
     rho0_NH3 = 0.7713  #g/L
     rho0_H2 = 0.0899  #g/L
     gas_properties = {
    'air': {'rho': rho0_air},
    'H2': {'rho': rho0_H2},
    'CH4': {'rho': rho0_CH4},
    'NH3': {'rho': rho0_NH3}
     }
     rho = gas_properties[species]['rho']
     debit_kgs = (Flow_rate_NL * rho) / (60*1000)
     return debit_kgs

def Kg_TO_NormoLitre(Flow_rate_Kg,species):
     # masse volumique à 0°C et 1 atm — source : Air Liquide
     rho0_air = 1.292  #g/L
     rho0_CH4 = 0.7173 #g/L
     rho0_NH3 = 0.7713  #g/L
     rho0_H2 = 0.0899  #g/L
     gas_properties = {
    'air': {'rho': rho0_air},
    'H2': {'rho': rho0_H2},
    'CH4': {'rho': rho0_CH4},
    'NH3': {'rho': rho0_NH3}
     }
     rho = gas_properties[species]['rho']
     debit_NLmin = (Flow_rate_Kg*60*1000) /rho
     return debit_NLmin

def Y_O2air():
     W_O2 = 0.032 #kg/mol
     W_N2 = 0.028 #kg/mol
     n_O2 = 0.21    #mol
     n_N2 = 1-n_O2  #mol
     Y_O2 = (n_O2 * W_O2)/(n_O2 * W_O2 + n_N2 * W_N2)
     return Y_O2

def stoich_rapport(fuel):
     gas = ct.Solution("gri30.yaml")
     s = gas.stoich_air_fuel_ratio(fuel=fuel,oxidizer="O2:1")
     return s

def massAIR_To_massFuel_in_premixed(m_AIR,Phi,fuel):
     s = stoich_rapport(fuel)
     Y_O2 = Y_O2air()
     m_F = Phi/s * m_AIR * Y_O2
     return m_F

def massF_FROM_Thermic_POWER(Power, PCI):
     Power = Power *1000
     PCI= PCI *1e6
     return Power/PCI

def CHECK_Globale_Equivalence_ratio(m_AIR,PhiGlobal,m_NH3,m_H2):
     W_H2 = 0.002 #kg/mol
     W_NH3 = 0.017 #kg/mol
     Beta = (m_H2/W_H2)/((m_H2/W_H2)+(m_NH3/W_NH3))
     fuel = f"NH3:{1-Beta},H2:{Beta}"
     sgl = stoich_rapport(fuel)
     print(f'Beta from the blend NH3 H2 is {Beta:3f}')
     print(f"\n The fuel composition is NH3: {1-Beta:3f}, H2: {Beta:3f}")
     print(f"The stochiometric rapport for this fuel is {sgl}")
     Y_O2 = Y_O2air()
     Phi = sgl * (m_H2 + m_NH3)/(m_AIR * Y_O2)
     print(f"The Equivalence ratio calculate is {Phi:3f}, the one form the exp is {PhiGlobal}")
     print()

def Mass_fraction_Injection_HYLON_Composition(m_AIR,m_H2,m_NH3,T,P,PhiClow):
     print("WARNING: This composition is valided with the follow HYLON configuration  ")
     print("\t WARNING: Hydrogen in the central jet ")
     print("\t WARNING: Ammonia mixing with air in the coflow ")

     W_H2  = 0.002
     W_NH3 = 0.017
     W_AIR = 0.21 * 0.032 + 0.79 * 0.028

     print("\n" + "="*45)
     print("          MASS FLOW RATES")
     print("="*45)
     print(f"{'Species':<8} | {'mdot [g/s]':>15}. | {'mdot [NL/min]':>15}")
     print("-"*45)
     print(f"{'H2':<8} | {m_H2*1000:15.3f} | {Kg_TO_NormoLitre(m_H2,'H2'):15.2f}")
     print(f"{'NH3':<8} | {m_NH3*1000:15.3f} | {Kg_TO_NormoLitre(m_NH3,'NH3'):15.2f}")
     print(f"{'Air':<8} | {m_AIR*1000:15.3f} | {Kg_TO_NormoLitre(m_AIR,'air'):15.2f}")
     print("="*45 + "\n")

     Y = {"H2": 1.0, "NH3": 0.0, "N2": 0.0, "O2": 0.0}
     print("\n--- MASS FRACTIONS Central jet Hydrogen ---")
     for k, v in Y.items():
          print(f"  • {k:<4} : {v:10.6f}")
          print("----------------------\n")

     Y_O2 = Y_O2air()
     Y_NH3_cf = m_NH3/(m_NH3+m_AIR)
     Y_O2_cf = (m_AIR * Y_O2)/(m_NH3+m_AIR)
     Y_N2_cf = (m_AIR * (1-Y_O2))/(m_NH3+m_AIR)

     Y = {"H2": 0.0, "NH3": Y_NH3_cf, "N2": Y_N2_cf, "O2": Y_O2_cf}
     print("\n--- MASS FRACTIONS CoFlow ---")
     for k, v in Y.items():
          print(f"  • {k:<4} : {v:10.6f}")
          print("----------------------\n")
     fuel = f"NH3:1.0"
     snh3 = stoich_rapport(fuel)
     Phi_Nh3 = snh3 * Y_NH3_cf/Y_O2_cf
     print(f"The Equivalence ratio in the Coflow calculate is {Phi_Nh3}, the one form the exp is {PhiClow}")
     

def __main__():
     PhiGlobal = 0.29
     PhiClow = 0.096
     alpha_Ph2 = 0.7
     Thermic_Power = 8 # KW
     PCI_H2 = 120 #MJ/kg
     PCI_NH3 = 18.6 #MJ/kg
     P = 101325  #Pa
     T = 273 #K
     m_AIR_NL = 378

     m_AIR_Kg = NormoLitre_TO_kg(m_AIR_NL,'air')

     m_NH3 = massAIR_To_massFuel_in_premixed(m_AIR_Kg,PhiClow,"NH3:1.0")

     m_H2 = massF_FROM_Thermic_POWER(Thermic_Power*alpha_Ph2, PCI_H2)

     CHECK_Globale_Equivalence_ratio(m_AIR_Kg,PhiGlobal,m_NH3,m_H2)

     Mass_fraction_Injection_HYLON_Composition(m_AIR_Kg,m_H2,m_NH3,T,P,PhiClow)

__main__()