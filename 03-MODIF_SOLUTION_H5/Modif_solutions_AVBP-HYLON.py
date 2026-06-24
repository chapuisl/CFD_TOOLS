"""
Requires : xmf_writer.py dans le même dossier
======================
Modifie les fractions massiques des espèces dans un fichier solution CFD HDF5
 en définissant un volume cylindrique (plein ou creux).

Fonctionnalités :
  - Lecture des coordonnées depuis le fichier mesh .h5
  - Sélection d'un volume cylindrique (axe X, Y ou Z)
  - Modification interactive des fractions massiques Y_k
  - Renormalisation automatique des espèces non modifiées
  - Mise à jour cohérente de rho, rhoE, rhou, rhov, rhow, r_bar
  - Gradient de transition à la frontière du volume (smoothstep)

Hypothèses solveur (YALES2) :
  - Variables conservatives : rho, rhoE, rhou, rhov, rhow
  - Espèces stockées en rhoY_k dans RhoSpecies/
  - r_bar = R_universal / M_mix  [J/kg/K]
  - rhoE = rho * e_tot,  e_tot = e_int + 0.5*|u|^2
  - e_int conservé lors du changement de composition (T fixée)
  - Loi des gaz parfaits : rho = P / (r_bar * T)
"""

import h5py
import numpy as np
import shutil
import os
import sys

# XMF writer (même dossier)
try:
    from xmf_writer import generate_xmf
    XMF_AVAILABLE = True
except ImportError:
    XMF_AVAILABLE = False

# ─────────────────────────────────────────────
#  CONSTANTE
# ─────────────────────────────────────────────
R_UNIVERSAL = 8.314  # J/mol/K

# Masses molaires [kg/mol] — à compléter si d'autres espèces présentes
MOLAR_MASSES = {
    "H2":   2.01588e-3,
    "N2":   28.0134e-3,
    "O2":   31.9988e-3,
    "NH3":  17.0305e-3,
    "H2O":  18.0153e-3,
    "CO":   28.0101e-3,
    "CO2":  44.0095e-3,
    "CH4":  16.0425e-3,
    "OH":   17.0079e-3,
    "NO":   30.0061e-3,
    "AR":   39.9480e-3,
    "HE":    4.0026e-3,
}


# ═══════════════════════════════════════════════════════════════════
#  1.  LECTURE DU MAILLAGE
# ═══════════════════════════════════════════════════════════════════

def load_mesh_coords(mesh_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Charge les coordonnées nodales X, Y, Z depuis le fichier mesh HDF5.

    Le chemin exact des coordonnées dépend de ton solveur/version.
    Le script tente plusieurs chemins classiques YALES2 et affiche
    la structure si aucun ne convient.
    """
    CANDIDATE_PATHS = [
        # (groupe, dataset_x, dataset_y, dataset_z)
        ("Coordinates",        "x",  "y",  "z"),
        ("Coordinates",        "X",  "Y",  "Z"),
        ("Nodes",              "x",  "y",  "z"),
        ("mesh/Coordinates",   "x",  "y",  "z"),
        ("Geometry/Nodes",     "x",  "y",  "z"),
        ("",                   "x",  "y",  "z"),   # à la racine
        ("",                   "X",  "Y",  "Z"),
    ]

    with h5py.File(mesh_file, "r") as f:
        for (group, dx, dy, dz) in CANDIDATE_PATHS:
            try:
                base = f[group] if group else f
                x = base[dx][:]
                y = base[dy][:]
                z = base[dz][:]
                print(f"  ✓ Coordonnées trouvées dans '{group or '/'}' : {dx}, {dy}, {dz}")
                print(f"    Nombre de nœuds : {len(x):,}")
                return x.astype(np.float64), y.astype(np.float64), z.astype(np.float64)
            except (KeyError, TypeError):
                continue

        # Aucun chemin trouvé → affiche la structure pour diagnostic
        print("\n⚠ Coordonnées non trouvées automatiquement.")
        print("  Structure du fichier mesh :")
        def _show(name, obj):
            tag = "DATASET" if isinstance(obj, h5py.Dataset) else "GROUP"
            shape = f" shape={obj.shape}" if isinstance(obj, h5py.Dataset) else ""
            print(f"    [{tag}] {name}{shape}")
        f.visititems(_show)
        raise KeyError(
            "Impossible de localiser les coordonnées XYZ.\n"
            "Renseigne manuellement CANDIDATE_PATHS dans load_mesh_coords()."
        )


# ═══════════════════════════════════════════════════════════════════
#  2.  SÉLECTION DES NŒUDS DANS LE VOLUME CYLINDRIQUE
# ═══════════════════════════════════════════════════════════════════

def radial_distance(x, y, z, center, axis):
    """Distance radiale par rapport à l'axe du cylindre."""
    cx, cy, cz = center
    if axis == "z":
        return np.sqrt((x - cx) ** 2 + (y - cy) ** 2), z - cz
    elif axis == "x":
        return np.sqrt((y - cy) ** 2 + (z - cz) ** 2), x - cx
    elif axis == "y":
        return np.sqrt((x - cx) ** 2 + (z - cz) ** 2), y - cy
    else:
        raise ValueError(f"Axe invalide : '{axis}'. Choisir parmi x, y, z.")


def select_nodes_cylinder(
    x, y, z,
    center: tuple,
    axis: str,
    r_inner: float,
    r_outer: float,
    height: float,
    h_offset: float = 0.0,
) -> np.ndarray:
    """
    Retourne les indices des nœuds dans le cylindre défini.

    Parameters
    ----------
    center   : (x0, y0, z0) — point de départ de l'axe (base du cylindre)
    axis     : 'x', 'y' ou 'z'
    r_inner  : rayon intérieur  (0 pour cylindre plein)
    r_outer  : rayon extérieur
    height   : longueur du cylindre le long de l'axe
    h_offset : décalage de la base le long de l'axe (0 par défaut)
    """
    r, along = radial_distance(x, y, z, center, axis)
    mask = (
        (r  >= r_inner)         &
        (r  <= r_outer)         &
        (along >= h_offset)     &
        (along <= h_offset + height)
    )
    idx = np.where(mask)[0]
    return idx


# ═══════════════════════════════════════════════════════════════════
#  3.  ZONE DE GRADIENT À LA FRONTIÈRE (SMOOTHSTEP)
# ═══════════════════════════════════════════════════════════════════

def smootherstep(t):
    t = np.clip(t,0,1)
    return t*t*t*(t*(t*6 - 15) + 10)


def build_transition_layers(
    x, y, z,
    center, axis,
    r_inner, r_outer, height,
    h_offset=0.0,
    n_layers: int = 8,
    transition_thickness: float = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    Construit les couches de transition à la frontière externe du cylindre.

    Retourne une liste de (indices_noeuds, poids_alpha) où :
      alpha = 0  → valeur imposée (intérieur)
      alpha = 1  → valeur originale conservée (extérieur)

    Parameters
    ----------
    n_layers             : nombre de couches de transition
    transition_thickness : épaisseur totale de la zone de transition.
                           Par défaut = r_outer * 0.10  (10 % du rayon)
    """
    if transition_thickness is None:
        transition_thickness = r_outer * 0.10

    r, along = radial_distance(x, y, z, center, axis)

    # Masque de base : zone légèrement élargie
    r_max_ext = r_outer + transition_thickness
    h_lo = h_offset - transition_thickness
    h_hi = h_offset + height + transition_thickness

    # Nœuds candidats à la transition (ni strictement intérieurs, ni loin)
    in_axial   = (along >= h_lo) & (along <= h_hi)
    in_radial  = (r <= r_max_ext)
    # Exclure les nœuds déjà dans le cœur du volume
    in_core    = (
        (r  >= r_inner) & (r  <= r_outer) &
        (along >= h_offset) & (along <= h_offset + height)
    )
    candidate = (in_axial &in_radial & (r >= r_inner) & ~in_core
)
    # Pour chaque nœud candidat, calcule la distance au bord le plus proche
    # Distance au bord radial extérieur (positive = dehors)
    d_r_outer = np.maximum(r - r_outer, 0.0)
    d_r_inner = np.maximum(r_inner - r, 0.0)
    # Distance au bord axial (positive = dehors)
    d_h_lo = np.maximum(h_offset - along, 0.0)             # positif si avant la base
    d_h_hi  = np.maximum(along - (h_offset + height), 0.0)    # positif si après le sommet
    d_axial = np.maximum(d_h_lo, d_h_hi)

    # Distance au bord = max des distances signées (>0 dehors)
    d_border = np.sqrt(d_r_outer**2 + d_r_inner**2 + d_axial**2)
    d_border = np.maximum(d_border, 0.0)
    idx_transition = np.where(candidate)[0]

    alpha_transition = np.clip(
        d_border[idx_transition] / transition_thickness,
        0.0,
        1.0
    )

    return idx_transition, alpha_transition
    # layers = []
    # for i in range(1, n_layers + 1):
    #     d_lo = (i - 1) * transition_thickness / n_layers
    #     d_hi = i       * transition_thickness / n_layers
    #     layer_mask = candidate & (d_border >= d_lo) & (d_border < d_hi)
    #     idx = np.where(layer_mask)[0]
    #     if len(idx) == 0:
    #         continue
    #     # Position normalisée dans la couche → alpha (0=intérieur, 1=extérieur)
    #     t = (d_border[idx] - d_lo) / (d_hi - d_lo)
    #     alpha = smootherstep(d_border / transition_thickness)
    #     layers.append((idx, alpha))

    # return layers


# ═══════════════════════════════════════════════════════════════════
#  4.  RECALCUL DES VARIABLES CONSERVATIVES
# ═══════════════════════════════════════════════════════════════════

def compute_r_bar(Y: dict, species_list: list) -> np.ndarray:
    """
    r_bar = R_univ * Σ(Y_k / M_k)   [J/kg/K]
    Constante spécifique du mélange.
    """
    inv_M_mix = np.zeros_like(list(Y.values())[0], dtype=np.float64)
    for sp in species_list:
        if sp not in MOLAR_MASSES:
            raise KeyError(
                f"Masse molaire inconnue pour '{sp}'. "
                f"Ajoute-la dans le dictionnaire MOLAR_MASSES."
            )
        inv_M_mix += Y[sp] / MOLAR_MASSES[sp]
    return R_UNIVERSAL * inv_M_mix


# ═══════════════════════════════════════════════════════════════════
#  5.  MODIFICATION PRINCIPALE Species
# ═══════════════════════════════════════════════════════════════════

def modify_species(
    output_file,
    idx_core,
    species_changes,
    idx_transition,
    alpha_transition,
):
    """
    Applique les nouvelles fractions massiques sur les nœuds idx_core,
    puis applique un gradient smoothstep sur les couches de transition.

    Parameters
    ----------
    output_file       : chemin du fichier HDF5 de sortie (copie modifiable)
    idx_core          : indices des nœuds dans le cœur du volume
    species_changes   : {nom_espece: nouvelle_valeur_Y}  (les autres sont renormalisées)
    transition_layers : liste de (idx_layer, alpha) produite par build_transition_layers
    """
    with h5py.File(output_file, "r+") as f:
        species_list = list(f["RhoSpecies"].keys())

        # ── Lire TOUTES les Y actuelles ──────────────────────────────
        rho_orig = f["GaseousPhase/rho"][:]
        Y_orig   = {sp: f[f"RhoSpecies/{sp}"][:] / rho_orig for sp in species_list}

        # ── Construire Y_new pour les nœuds du cœur ─────────────────
        Y_new = {sp: Y_orig[sp].copy() for sp in species_list}

        fixed_species = list(species_changes.keys())
        fixed_sum     = sum(species_changes.values())

        if fixed_sum > 1.0 + 1e-10:
            raise ValueError(
                f"Somme des fractions massiques imposées = {fixed_sum:.4f} > 1 !"
            )

        free_species     = [sp for sp in species_list if sp not in fixed_species]
        free_sum_at_core = sum(Y_orig[sp][idx_core] for sp in free_species)
        remaining        = 1.0 - fixed_sum

        # Appliquer les nouvelles Y imposées
        for sp in fixed_species:
            Y_new[sp][idx_core] = species_changes[sp]

        # Renormaliser les espèces libres
        for sp in free_species:
            denom = free_sum_at_core
            if np.any(denom > 0):
                ratio = np.where(denom > 0,
                                 Y_orig[sp][idx_core] / denom * remaining,
                                 remaining / len(free_species))
            else:
                ratio = np.full(len(idx_core), remaining / len(free_species))
            Y_new[sp][idx_core] = ratio

        # ── Appliquer le gradient sur les couches de transition ──────
        # alpha=0 → valeur du cœur modifié, alpha=1 → valeur originale
        #
        # Pour chaque couche, la "valeur cœur" est interpolée à partir
        # de la valeur imposée (species_changes) et des originaux renormalisés.
        # On simplifie : la valeur cible pour alpha=0 est celle du cœur.

        # Valeur imposée constante côté intérieur
        Y_inner = {sp: species_changes.get(sp, None) for sp in species_list}

        # Calculer Y_inner pour les espèces libres (renormalisées globalement)
        free_sum_global = {
            sp: float(Y_orig[sp][idx_core].mean()) for sp in free_species
        }
        total_free_global = sum(free_sum_global.values())
        for sp in free_species:
            if total_free_global > 0:
                Y_inner[sp] = free_sum_global[sp] / total_free_global * remaining
            else:
                Y_inner[sp] = remaining / len(free_species)

        for sp in species_list:
            Y_new[sp][idx_transition] = (
                (1.0 - alpha_transition) * Y_inner[sp]
                + alpha_transition * Y_orig[sp][idx_transition]
            )

        # ── Combiner cœur + couches de transition ────────────────────
        all_modified_idx = np.unique(
    np.concatenate([
        idx_core,
        idx_transition
    ])
)
        # ── Mise à jour des variables conservatives ──────────────────
        rho  = f["GaseousPhase/rho"][:]
        rhoE = f["GaseousPhase/rhoE"][:]
        rhou = f["GaseousPhase/rhou"][:]
        rhov = f["GaseousPhase/rhov"][:]
        rhow = f["GaseousPhase/rhow"][:]
        P    = f["Additionals/pressure"][:].astype(np.float64)
        T    = f["Additionals/temperature"][:].astype(np.float64)

        idx_m = all_modified_idx

        # Énergie totale massique AVANT modification (conservée)
        e_tot_old = rhoE[idx_m] / rho[idx_m]

        # Vitesse AVANT modification (conservée)
        u_old = rhou[idx_m] / rho[idx_m]
        v_old = rhov[idx_m] / rho[idx_m]
        w_old = rhow[idx_m] / rho[idx_m]

        # Nouveau r_bar et rho
        Y_idx      = {sp: Y_new[sp][idx_m] for sp in species_list}
        r_bar_new  = compute_r_bar(Y_idx, species_list)
        rho_new    = P[idx_m] / (r_bar_new * T[idx_m])

        # Écriture des nouvelles rhoY
        for sp in species_list:
            rhoY = f[f"RhoSpecies/{sp}"][:]
            rhoY[idx_m] = Y_new[sp][idx_m] * rho_new
            f[f"RhoSpecies/{sp}"][:] = rhoY

        # Écriture de rho
        rho[idx_m] = rho_new
        f["GaseousPhase/rho"][:] = rho

        # rhoE : e_tot massique conservé
        rhoE[idx_m] = rho_new * e_tot_old
        f["GaseousPhase/rhoE"][:] = rhoE

        # rhou/rhov/rhow : vitesse conservée
        rhou[idx_m] = rho_new * u_old
        rhov[idx_m] = rho_new * v_old
        rhow[idx_m] = rho_new * w_old
        f["GaseousPhase/rhou"][:] = rhou
        f["GaseousPhase/rhov"][:] = rhov
        f["GaseousPhase/rhow"][:] = rhow

        # Mettre à jour r_bar dans Additionals (si présent)
        if "r_bar" in f["Additionals"]:
            r_bar_arr = f["Additionals/r_bar"][:]
            r_bar_arr[idx_m] = r_bar_new.astype(np.float32)
            f["Additionals/r_bar"][:] = r_bar_arr

        # Résumé
        n_core  = len(idx_core)
        n_trans = len(all_modified_idx) - n_core
        print(f"\n  ✓ Nœuds cœur modifiés       : {n_core:,}")
        print(f"  ✓ Nœuds transition modifiés  : {n_trans:,}")
        print(f"  ✓ Total nœuds mis à jour     : {len(all_modified_idx):,}")
        print(f"  ✓ rho_new  : [{rho_new.min():.4f} – {rho_new.max():.4f}] kg/m³")


# ═══════════════════════════════════════════════════════════════════
#  6.  VÉRIFICATION POST-MODIFICATION
# ═══════════════════════════════════════════════════════════════════

def verify_solution(output_file: str, idx_core: np.ndarray):
    """Affiche un bilan de vérification après modification."""
    print("\n─── Vérification de la solution modifiée ───")
    with h5py.File(output_file, "r") as f:
        rho      = f["GaseousPhase/rho"][:]
        species  = list(f["RhoSpecies"].keys())

        # Somme des fractions massiques
        sum_Y = sum(f[f"RhoSpecies/{sp}"][:] for sp in species) / rho
        print(f"  ΣY_k  — global  : [{sum_Y.min():.6f} – {sum_Y.max():.6f}]  (doit être ~1)")
        print(f"  ΣY_k  — cœur    : [{sum_Y[idx_core].min():.6f} – {sum_Y[idx_core].max():.6f}]")

        print(f"\n  rho   — global  : [{rho.min():.4f} – {rho.max():.4f}] kg/m³")
        print(f"  rho   — cœur    : [{rho[idx_core].min():.4f} – {rho[idx_core].max():.4f}] kg/m³")

        print("\n  Fractions massiques dans le cœur :")
        for sp in species:
            Y_core = f[f"RhoSpecies/{sp}"][:][idx_core] / rho[idx_core]
            print(f"    Y_{sp:<6s}: [{Y_core.min():.6f} – {Y_core.max():.6f}]  "
                  f"moy={Y_core.mean():.6f}")


# ═══════════════════════════════════════════════════════════════════
#  7.  INTERFACE UTILISATEUR INTERACTIVE
# ═══════════════════════════════════════════════════════════════════

def prompt_float(msg, default=None):
    while True:
        raw = input(msg).strip()
        if raw == "" and default is not None:
            return default
        try:
            return float(raw)
        except ValueError:
            print("  ⚠ Entier ou décimal attendu. Réessaie.")


def prompt_choice(msg, choices):
    choices_lower = [c.lower() for c in choices]
    while True:
        raw = input(msg).strip().lower()
        if raw in choices_lower:
            return raw
        print(f"  ⚠ Choix valides : {choices}")


def main():
    print("=" * 60)
    print("   Modification de solution CFD — Espèces & Volume")
    print("=" * 60)

    # ── Fichiers d'entrée/sortie ─────────────────────────────────
    sol_file  = input("\nFichier solution   (.h5) : ").strip()
    mesh_file = input("Fichier maillage   (.h5) : ").strip()
    out_file  = input("Fichier de sortie  (.h5) [sol_modified.h5] : ").strip()
    if not out_file:
        out_file = "sol_modified.h5"

    if not os.path.isfile(sol_file):
        sys.exit(f"✗ Fichier solution introuvable : {sol_file}")
    if not os.path.isfile(mesh_file):
        sys.exit(f"✗ Fichier maillage introuvable : {mesh_file}")

    print(f"\nCopie de {sol_file} → {out_file} …")
    shutil.copy2(sol_file, out_file)

    # ── Chargement des coordonnées ───────────────────────────────
    print("\nChargement des coordonnées …")
    x, y, z = load_mesh_coords(mesh_file)
    n_nodes  = len(x)
    print(f"  Domaine X : [{x.min():.4f} – {x.max():.4f}]")
    print(f"  Domaine Y : [{y.min():.4f} – {y.max():.4f}]")
    print(f"  Domaine Z : [{z.min():.4f} – {z.max():.4f}]")

    # ── Paramètres du volume cylindrique ────────────────────────
    print("\n── Définition du volume cylindrique ──")
    print("  Centre de la BASE du cylindre (point 1) :")
    cx = prompt_float("    x0 : ")
    cy = prompt_float("    y0 : ")
    cz = prompt_float("    z0 : ")
    center = (cx, cy, cz)

    axis   = prompt_choice("  Axe du cylindre (x / y / z) : ", ["x", "y", "z"])
    height = prompt_float("  Hauteur du cylindre (le long de l'axe) : ")

    cyl_type = prompt_choice("  Type de cylindre (plein / creux) : ", ["plein", "creux"])
    r_outer  = prompt_float("  Rayon extérieur : ")
    r_inner  = 0.0
    if cyl_type == "creux":
        r_inner = prompt_float("  Rayon intérieur : ")
        if r_inner >= r_outer:
            sys.exit("✗ r_inner doit être < r_outer.")

    # ── Zone de gradient ─────────────────────────────────────────
    print("\n── Zone de gradient à la frontière ──")
    default_thick = 0.0005
    thick = prompt_float(
        f"  Épaisseur de la zone de transition [défaut={default_thick}] : ",
        default=default_thick
    )
    n_lay = int(prompt_float("  Nombre de couches de transition [défaut=8] : ", default=8))

    # ── Sélection des nœuds ──────────────────────────────────────
    print("\nRecherche des nœuds dans le volume …")
    idx_core = select_nodes_cylinder(x, y, z, center, axis, r_inner, r_outer, height)
    print(f"  → {len(idx_core):,} nœuds sélectionnés dans le cœur")
    if len(idx_core) == 0:
        sys.exit("✗ Aucun nœud trouvé. Vérifie les coordonnées et l'unité.")

    print("Construction des couches de transition …")
    #TODO
    # trans_layers = build_transition_layers(
    #     x, y, z, center, axis,
    #     r_inner, r_outer, height,
    #     n_layers=n_lay,
    #     transition_thickness=thick,
    # )

    idx_transition, alpha_transition = build_transition_layers(
        x, y, z, center, axis,
        r_inner, r_outer, height,
        n_layers=n_lay,
        transition_thickness=thick,
    )

    # n_trans_nodes = sum(len(i) for i, _ in trans_layers)
    # print(f"  → {n_trans_nodes:,} nœuds dans la zone de transition ({n_lay} couches)")

    # ── Choix des espèces à modifier ────────────────────────────
    with h5py.File(out_file, "r") as f:
        species_list = list(f["RhoSpecies"].keys())

    print("\n── Espèces disponibles ──")
    for i, sp in enumerate(species_list, 1):
        print(f"  {i}. {sp}")

    species_changes = {}
    print("\nSaisis les espèces à modifier (Entrée sans valeur = terminer) :")
    while True:
        sp_input = input("  Espèce (nom exact) : ").strip()
        if sp_input == "":
            break
        if sp_input not in species_list:
            print(f"  ⚠ '{sp_input}' non trouvé. Espèces valides : {species_list}")
            continue
        val = prompt_float(f"  Y_{sp_input} = ")
        if val < 0.0 or val > 1.0:
            print("  ⚠ Fraction massique doit être entre 0 et 1.")
            continue
        species_changes[sp_input] = val

    if not species_changes:
        sys.exit("✗ Aucune espèce sélectionnée. Abandon.")

    # Vérification somme
    total = sum(species_changes.values())
    print(f"\n  Somme des Y imposées : {total:.4f}")
    if total > 1.0:
        sys.exit(f"✗ Somme > 1.0  ({total:.4f}). Corrige les valeurs.")
    remaining = 1.0 - total
    free = [sp for sp in species_list if sp not in species_changes]
    print(f"  Reste à distribuer sur {free} : {remaining:.4f}")

    # ── Récapitulatif avant modification ─────────────────────────
    print("\n" + "─" * 50)
    print("  RÉCAPITULATIF")
    print(f"  Fichier sortie   : {out_file}")
    print(f"  Centre           : {center}")
    print(f"  Axe / Hauteur    : {axis} / {height}")
    print(f"  Rayons           : r_inner={r_inner}, r_outer={r_outer}")
    print(f"  Nœuds cœur      : {len(idx_core):,}")
    print(f"  Transition       : épaisseur={thick}, couches={n_lay}")
    print(f"  Modifications    : {species_changes}")
    print("─" * 50)
    confirm = input("  Confirmer la modification ? (oui/non) : ").strip().lower()
    if confirm not in ("oui", "o", "yes", "y"):
        sys.exit("Annulé.")

    # ── Application 
    print("\nApplication des modifications …")
    modify_species(out_file, idx_core, species_changes, idx_transition, alpha_transition)

    # ── Vérification 
    verify_solution(out_file, idx_core)

    print(f"\n Fichier HDF5 sauvegardé : {out_file}")

    # ── Génération du XMF pour Paraview 
    if not XMF_AVAILABLE:
        print("\n⚠ xmf_writer.py non trouvé — XMF non généré.")
        print("  Place xmf_writer.py dans le même dossier que ce script.")
    else:
        print("\n── Export XMF pour Paraview ──")
        xmf_mode = prompt_choice(
            "  Mode XMF — dual (mesh + sol séparés) / single (tout dans sol) : ",
            ["dual", "single"]
        )
        # Temps physique depuis le fichier solution
        with h5py.File(out_file, "r") as f:
            t_phys = float(f["Parameters/dtsum"][0])

        xmf_out = os.path.splitext(out_file)[0] + ".xmf"
        generate_xmf(
            sol_file   = out_file,
            mesh_file  = mesh_file,
            xmf_file   = xmf_out,
            time_value = t_phys,
            mode       = xmf_mode,
        )
        print(f"\n Open '{xmf_out}' in Paraview for visualisation.")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    main()