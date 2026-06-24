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
from Global_plot_figure2 import plot_evolution
import Graphic_Configuration as GC
"""
# ===================================================================================================================
#  Import library
# ===================================================================================================================
"""
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib as mpl
import sys
import csv
import os 


"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def config_plot(style="default", title_size=10, label_size=30, tick_size=35, background="#f9f9f9"):
    plt.style.use(style)
    mpl.rcParams.update({
        "figure.facecolor": background,
        "axes.facecolor": "white",
        "axes.edgecolor": "#333333",
        "axes.labelsize": label_size,
        "axes.titlesize": title_size,
        "axes.titleweight": "bold",
        "grid.alpha": 1.0,
        "grid.color": "#888888",
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        # "legend.fontsize": 20,
        # "legend.frameon": True,
        # "legend.facecolor": "#f0f0f0",
        # "legend.edgecolor": "#cccccc",
        # "lines.linewidth": 4,
        # "axes.linewidth": 2.5,
        # "grid.linewidth": 2.0,
        # "lines.markersize": 25,
        # "lines.markeredgewidth": 2,
        "xtick.major.width": 6,
        "ytick.major.width": 6,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.pad": 6,
        "ytick.major.pad": 6,
        "xtick.minor.width": 1.5,
        "ytick.minor.width": 1.5,
        "xtick.minor.size": 5,  # Ajustez la taille des ticks mineurs ici
        "ytick.minor.size": 5, 
    })

def PLOT_INITIALISATION():
    COLOR = GC.Black_Blue_short()
    STYLE = GC.LineStyles()
    config_plot()
    return STYLE, COLOR

STYLE, COLOR = PLOT_INITIALISATION()

def plot(data_final,y_max = None):
    fig, ax_base = plt.subplots(figsize=(12, 7))
    axes = [ax_base]
    lines = []

    for i, d in enumerate(data_final):
        if i == 0:
            ax = ax_base
        else:
            ax = ax_base.twinx()
            ax.spines["right"].set_position(("outward", 130 * (i - 1)))  # décale chaque axe
        
        l, = ax.plot(Position, d["values"], color=COLOR[i],linewidth= 3.0)
        ax.set_ylabel(f'{d["label"]} ({d["unit"]})', color=COLOR[i])
        ax.tick_params(axis="y", labelcolor=COLOR[i])
        if y_max is not None:
            ax.set_ylim(0, y_max)  # ← appel correct, pas une assignation

        if i ==0:
            plt.grid(True)
        axes.append(ax)
        lines.append(l)

    # Légende commune sur le premier axe
    ax_base.set_xlabel("Time (s)")
    plt.tight_layout()
    plt.show()

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
    read_dat = read_dat_probe(probe)
    all_data = all_data | read_dat
    

# --- Plots : un par paramètre, une courbe par fichier ---

params = [col for col in all_data.keys() if col != 'x']
base_name = (
    f"{{}}_EVOLUTION_SPATIAL"
)

data_finalNS = []
data_final2NS = []

data_finalspec = []
data_final2spec = []
it = 0 

unity = {'T': ' [K]', 'P': ' [bar]','rho':' [kg/m3]','u':' [cm/s]', 'v':' [m/s]'}
ValNS = ['T','P','u']
Valspec = ['NH3','N2','O2','H2O']

for param in params:
    if param not in ValNS:
        continue
    Position = []
    data     = []
    labels   = []

    n=1
    Position = all_data['x'].tolist()
    data = all_data[param].tolist()
    
    if param == 'P':
        data = [val * 1e-5 for val in data] 
    
    if param == 'u':
        data = [val * 1e2 for val in data] 

    prefix = 'Y_'
    try:
        unit = unity[param]
        prefix = '' 
    except:
        unit = ''

    if param == 'T':
        dataT = [{"label": param, "values": data , "unit": unit}]
    else:
        data_finalNS.append({"label": param, "values": data , "unit": unit})
        
    it +=1
data_finalNS = dataT + data_finalNS

it = 0
for param in params:
    if param not in Valspec:
        continue
    Position = []
    data     = []
    labels   = []

    n=1
    Position = all_data['x'].tolist()
    data = all_data[param].tolist()
    
    prefix = 'Y_'
    try:
        unit = unity[param]
        prefix = '' 
    except:
        unit = ''

    if param == 'NH3':
        dataNH3 = [{"label": param, "values": data , "unit": unit}]
    else:
        data_finalspec.append({"label": param, "values": data , "unit": unit})
        
    it +=1
data_finalspec = dataNH3 + data_finalspec

plot(data_finalNS)
plot(data_finalspec,y_max = 1.0)

