"""
# ===================================================================================================================
#  Reading file
# ===================================================================================================================
#
#  Author          : Lilian CHAPUIS
#  Affiliation     : IMFT - Institut de Mécanique des Fluides de Toulouse
#  Location        : Toulouse, France
#  Creation Date   : 16  March 2026
#  Last Modified   : 16  March 2026
#  Version         : 1.0.01
#
# -------------------------------------------------------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------------------------------------------------------
#  Features:
#- **Converts AVBP temporal output files** to human-readable CSV format.
#- **Plots time evolution** of key statistics (min, max, mean) for selected variables, replicating the functionality of AVBP’s `xt` and `xm` tools.
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
    COLOR = GC.Black_Blue_short()
    STYLE = GC.LineStyles()
    GC.config_plot()
    return STYLE, COLOR

STYLE, COLOR = PLOT_INITIALISATION()

def read_dat_probe(probe_file):
    value = []
    header = []
    Start_bool = True
    NB_Same_line = 0
    with open(probe_file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in spamreader:
            # Filtrer les champs vides
            row_clean = list(filter(lambda a: a != '', row))
            if len(row_clean) == 0:
                continue
            if row_clean[0] == '#':
                # La ligne header est : # x ... ... 
                header = row_clean[1:]  # on ignore le '#'
            else:
                try:
                    val_save = [float(k) for k in row_clean]
                    if Start_bool is True:
                        value.append(val_save)
                        Start_bool = False
                    else:
                        if  val_save[0] == value[-1][0]:
                            NB_Same_line += 1
                            continue
                        else:
                            value.append(val_save)
                except ValueError:
                    continue  # ignorer les lignes non numériques
    print(f'{NB_Same_line} line was identical')
    val = np.transpose(np.array(value))
    data = dict()
    if len(header) == 0:
        header = ['col%i' % i for i in range(val.shape[0])]
    for i_h, head in enumerate(header):
        print(head)
        data[head] = val[i_h]

    return data


probe_list = sys.argv[1:]  # tous les fichiers passés en argument
file_name = sys.argv[1]
path = os.getcwd()

output_dir = os.path.join(path, "SOLUT_RESULTAT")
os.makedirs(output_dir, exist_ok=True)  

if len(probe_list) == 0:
    probe_list = sorted(Path('.').glob('full_probe*_readbin'))
    probe_list = [str(probe.name) for probe in probe_list]

all_data = {}
for probe in probe_list:
    all_data[probe] = read_dat_probe(probe)
    
# --- Sauvegarde CSV ---
for probe, probe_dict in all_data.items():
    probe_name = probe.split("/")[-1]

    data = pd.DataFrame.from_dict(probe_dict)
    data.to_csv(os.path.join(path, f'{probe_name}.csv'), index=False)
    print(f'PROBE DATA {probe} turned into csv.')

# --- Plots : un par paramètre, une courbe par fichier ---

sample = next(iter(all_data.values()))
params = [col for col in sample.keys() if col != 'x']

base_name = (
    f"{{}}_EVOLUTION_SPATIAL"
)
unity = {'T': ' [K]', 'P': ' [bar]','rho':' [kg/m3]','u':' [m/s]', 'v':' [m/s]'}
it = 0
for param in params:
    it += 1
    print(f'PLOT number {it}')
    Position = []
    data     = []
    labels   = []

    n=1
    for probe, probe_dict in all_data.items():
        if 'x' in probe_dict and param in probe_dict:
            Position.append(probe_dict['x'].tolist())
            data.append(probe_dict[param].tolist())
            labels.append(sys.argv[n])
            n+=1

    if param == 'P':
        data = [[val * 1e-5 for val in sublist] for sublist in data]

    if param == 'T':
        print('Temperature')
        for i in range(len(Position)):
            dTdx = np.gradient(data[i], Position[i])
            thickness = (np.max(data[i]) - np.min(data[i])) / np.max(np.abs(dTdx))
            print(f'grad max of dTdx {np.max(np.abs(dTdx))} for plot number {i} ')
            print(f'Thickness {thickness*1e3} mm for plot number {i} ')

    prefix = 'Y_'
    try:
        unit = unity[param]
        prefix = '' 
    except:
        unit = ''

    start = True
    if start is True:
        start = False
        plot_evolution(
        Position,
        data,
        labels,
        colors = COLOR,
        ylabel= prefix + param +unit,
        xlabel="x [m]",
        save_fig=True,
        save_path=output_dir,
        name_fig='00_start_label',
        # marker='o'
    )

    plot_evolution(
        Position,
        data,
        # labels,
        colors = COLOR,
        ylabel=prefix + param +unit,
        xlabel="x [m]",
        save_fig=True,
        save_path=output_dir,
        name_fig=base_name.format(param),
        # marker='o'
    )
          

    
