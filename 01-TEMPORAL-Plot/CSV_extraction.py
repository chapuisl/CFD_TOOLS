"""
# ===================================================================================================================
#  Reading file
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mécanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 9  March 2026
#  Last Modified   : 9  March 2026
#  Version         : 1.0.01
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  This python transform les fichier de sortie des temporaux de AVBP en CSV lisible.
# En plus de ca il plot les evolution min max et mean des valeurs donné dans les fichier temporaux 
# comme l'outil xt et xm de AVBP.
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
#  Imports from other modules
# ===================================================================================================================
"""
from Global_plot_figure import plot_evolution
import Graphic_Configuration as GC
"""
# ===================================================================================================================
#  Import library
# ===================================================================================================================
"""

from pathlib import Path
import pandas as pd
import numpy as np
import sys
import csv
import os 


"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""
def PLOT_INITIALISATION():
    COLOR = [(0.0, 0.2, 0.9), 
            (0.7, 0.1, 0.0),                                               
            (0.0, 0.0, 0.0),  ]     
    STYLE = GC.LineStyles()
    GC.config_plot()
    return STYLE, COLOR
     
STYLE, COLOR = PLOT_INITIALISATION()

def read_dat_probe(probe_file):
    value = []
    header = []
    with open(probe_file,newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in spamreader:
            if row[0] == r'#': header.append(row[1])
            else: value.append([float(k) for k in list(filter(lambda a: a != '', row)) ])
    val = np.transpose(np.array(value))
    data = dict()
    if len(header) == 0: header = ['col%i'%(i) for i in range(val.shape[0])]
    for i_h,head in enumerate(header): data[head] = val[i_h]

    return data


path = os.getcwd()
output_dir = os.path.join(path, "TEMPORAL_RESULTAT")
os.makedirs(output_dir, exist_ok=True)  

try:
    file_name = sys.argv[1]
    probe_list = [file_name]
except IndexError:
    probe_list = sorted(Path('.').glob('full_probe*_readbin'))
    probe_list = [str(probe.name) for probe in probe_list]
    
Species_list =[]
n=0
for probe in probe_list:
    
    probe_dict = read_dat_probe(probe)
    for name, val in probe_dict.items():
        n+=1
        print(f"Field {n}:",name)
        if name.startswith('Y') and name.endswith('min'):
            name_species = name.split("_")[1]
            Species_list.append(name_species)
            
        probe_dict[name] = val.tolist()

    data = pd.DataFrame.from_dict(probe_dict)
    data.to_csv(os.path.join(path, f'{probe}.csv'), index=False)
    
    print(f'PROBE DATA {probe} turned into csv.')
    
    
Time = probe_dict['atime']

Pressure_Evol = [
    [x * 1e-5 for x in probe_dict['P_mean']],
    [x * 1e-5 for x in probe_dict['P_max']],
    [x * 1e-5 for x in probe_dict['P_min']]
]

Temperature_Evol = [probe_dict['T_mean'],probe_dict['T_max'],probe_dict['T_min']]

Density_Evol = [probe_dict['rho_mean'],probe_dict['rho_max'],probe_dict['rho_min']]

HRR_Evol = [
    [x * 1e-9 for x in probe_dict['HR_mean']],
    [x * 1e-9 for x in probe_dict['HR_max']],
    [x * 1e-9 for x in probe_dict['HR_min']]
]

Norm_vel_Evol = [probe_dict['norm_vel_mean'],probe_dict['norm_vel_max'],probe_dict['norm_vel_min']]


print("\n")
print(" -->  SAVE PLOT")

base_name = (
    f"{{}}_EVOLUTION_TEMPORAL"
)
extension_label = ['_mean','_max','_min']

plots = [
    (Pressure_Evol, 'P', "P [bar]", "Pressure"),
    (Temperature_Evol, 'T', "T [K] ", "Temperature"),
    (Density_Evol, 'Rho', "rho [kg/m3] ", "Density"),
    (HRR_Evol, 'HRR', "HRR [GW/m3]", "Heat_Release"),
    (Norm_vel_Evol, 'Norm_vel', "Norm velo [m/s]", "Norm_vel"),
]

for spec in Species_list:
    data  = []
    for ext in extension_label:
        data.append(probe_dict['Y_'+ spec + ext])
    plots.append((data,spec,f"Y_ {spec}",f"Y_{spec}" ))
        
start = True
for data, extension, ylabel, suffix in plots:
    if start is True:
        start = False
        labels = [extension+ b for b in extension_label ]
        plot_evolution(
        Time,
        data,
        labels,
        colors = COLOR,
        ylabel=ylabel,
        xlabel="Time [s]",
        save_fig=True,
        save_path=output_dir,
        name_fig='00_start_label' + base_name.format(suffix),
    )
    plot_evolution(
        Time,
        data,
        # labels,
        colors = COLOR,
        ylabel=ylabel,
        xlabel="Time [s]",
        save_fig=True,
        save_path=output_dir,
        name_fig=base_name.format(suffix),
    )
print(" -->  END SAVE PLOT")















