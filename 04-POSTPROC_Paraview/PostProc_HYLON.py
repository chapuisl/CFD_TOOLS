"""
===============================================================================
 Script de post-traitement automatique ParaView - Bruleur HYLON (ABVP)
===============================================================================

Objectif :
    Automatiser tout ce que tu refais manuellement à chaque session ParaView :
        - fond blanc (pas de dégradé gris/bleu)
        - colormaps cohérentes par champ physique (T, vitesse, pression, Y_i...)
        - bornes de colorbar propres (fixes ou auto sur les données)
        - police et couleurs de la scalar bar lisibles sur fond blanc
        - orientation/format de la légende
        - export automatique de screenshots pour chaque champ

Utilisation :
    1. Ouvre ParaView.
    2. Charge ton fichier (.foam, .vtu, .pvd, etc.) normalement à la souris,
       OU adapte la fonction load_source() ci-dessous pour le faire charger
       automatiquement depuis ce script.
    3. View > Python Shell (ou Tools > Python Shell), puis :
           exec(open("posttraitement_hylon.py").read())
       OU en ligne de commande :
           pvpython posttraitement_hylon.py --input mon_cas.foam --outdir ./figures

    Le script détecte la source active si elle existe déjà dans ParaView,
    sinon il charge le fichier donné en argument.

-------------------------------------------------------------------------------
 A ADAPTER À TON CAS
-------------------------------------------------------------------------------
Modifie uniquement la section "CONFIGURATION" plus bas :
    - FIELDS : liste des champs à traiter, avec leur colormap et leurs bornes
    - CAMERA : vue par défaut (tu peux aussi la laisser telle quelle)
    - EXPORT : résolution et dossier de sortie des images

===============================================================================
"""

import os
import sys
import argparse

from paraview.simple import *
from paraview import servermanager as sm

paraview.simple._DisableFirstRenderCameraReset()


# =============================================================================
# CONFIGURATION - à adapter à ton cas HYLON / ABVP
# =============================================================================

# Dictionnaire des champs à mettre en forme.
# Pour chaque champ : nom exact tel qu'il apparaît dans ParaView,
#                      preset de colormap ParaView,
#                      bornes (min, max) ou None pour auto-rescale sur les données,
#                      échelle log (True/False),
#                      nombre de décimales affichées dans la légende.
FIELDS = {
    "T":        {"preset": "Cool to Warm",              "range": (280, 2400), "log": False, "decimals": 0, "label": "Temperature [K]"},
    "U":        {"preset": "Viridis (matplotlib)",       "range": None,        "log": False, "decimals": 1, "label": "Velocity magnitude [m/s]"},
    "p":        {"preset": "Cool to Warm",              "range": None,        "log": False, "decimals": 0, "label": "Pressure [Pa]"},
    "CH4":      {"preset": "Viridis (matplotlib)",       "range": (0, 1),      "log": False, "decimals": 2, "label": "Y_CH4 [-]"},
    "O2":       {"preset": "Viridis (matplotlib)",       "range": (0, 0.23),   "log": False, "decimals": 2, "label": "Y_O2 [-]"},
    "CO2":      {"preset": "Viridis (matplotlib)",       "range": (0, 0.2),    "log": False, "decimals": 2, "label": "Y_CO2 [-]"},
    "rho":      {"preset": "X Ray",                      "range": None,        "log": False, "decimals": 2, "label": "Density [kg/m3]"},
    "OH":       {"preset": "Inferno (matplotlib)",       "range": None,        "log": False, "decimals": 4, "label": "Y_OH [-]"},
    "vorticity":{"preset": "Cool to Warm",               "range": None,        "log": False, "decimals": 0, "label": "Vorticity [1/s]"},
}

# Couleur de fond de la vue (blanc pur)
BACKGROUND_RGB = [1.0, 1.0, 1.0]

# Apparence de la légende (scalar bar) - lisible sur fond blanc
LEGEND_FONT_COLOR = [0.0, 0.0, 0.0]   # noir
LEGEND_FONT_SIZE = 14
LEGEND_TITLE_FONT_SIZE = 16
LEGEND_ORIENTATION = "Horizontal"     # "Horizontal" ou "Vertical"
LEGEND_POSITION = [0.30, 0.03]        # position normalisée (x, y) bas de l'écran

# Options générales de rendu
SHOW_EDGES = False
SHOW_ORIENTATION_AXES = False
SHOW_CUBE_AXES = False

# Résolution des images exportées
EXPORT_RESOLUTION = [1920, 1080]

# Dossier de sortie par défaut
DEFAULT_OUTDIR = "./figures_hylon"


# =============================================================================
# FONCTIONS UTILITAIRES - normalement pas besoin d'y toucher
# =============================================================================

def get_active_source_or_load(input_path=None):
    """Récupère la source déjà chargée dans ParaView, sinon charge input_path."""
    src = GetActiveSource()
    if src is not None:
        print(f"[INFO] Source active détectée : {src}")
        return src

    if input_path is None:
        raise RuntimeError(
            "Aucune source active dans ParaView et aucun --input fourni. "
            "Charge un fichier dans ParaView avant d'exécuter ce script, "
            "ou lance avec pvpython --input chemin/vers/fichier"
        )

    ext = os.path.splitext(input_path)[1].lower()
    print(f"[INFO] Chargement du fichier : {input_path}")

    if ext == ".foam":
        src = OpenFOAMReader(FileName=input_path)
    elif ext in (".xmf", ".xdmf"):
        # ABVP écrit les chemins HDF5 (mesh, solution) en RELATIF par rapport
        # à l'emplacement du fichier .xmf lui-même. Le lecteur XDMF legacy de
        # ParaView, lui, résout ces chemins relatifs par rapport au cwd du
        # process -> si le script n'est pas lancé depuis le dossier exact du
        # .xmf, la résolution part n'importe où (cf. erreurs HDF5 "Unable to
        # open file"). On reproduit donc le comportement de l'IHM ParaView en
        # se plaçant temporairement dans le dossier du .xmf avant ouverture.
        input_abs = os.path.abspath(input_path)
        input_dir = os.path.dirname(input_abs)
        input_name = os.path.basename(input_abs)
        old_cwd = os.getcwd()
        try:
            os.chdir(input_dir)
            print(f"[INFO] cwd temporaire pour résolution des chemins relatifs XDMF : {input_dir}")
            src = XDMFReader(FileNames=[input_name])
            src.UpdatePipeline()
        finally:
            os.chdir(old_cwd)
        return src
    elif ext in (".vtu", ".vtk", ".vtm", ".vtp"):
        src = OpenDataFile(input_path)
    elif ext == ".pvd":
        src = PVDReader(FileName=input_path)
    else:
        src = OpenDataFile(input_path)

    src.UpdatePipeline()
    return src


def goto_timestep(source, timestep_index, xdmf_dir=None):
    """Positionne l'animation sur un pas de temps donné.
    timestep_index = -1 -> dernier instant disponible (défaut, typique pour un champ converge)
    timestep_index = 0  -> premier instant
    timestep_index = N  -> index N dans la liste des TimestepValues

    xdmf_dir : si fourni, on se replace temporairement dans ce dossier pendant
               le changement de temps, car changer AnimationTime redéclenche
               la lecture du .h5 avec les mêmes chemins relatifs que
               get_active_source_or_load() -> même contrainte de cwd.
    """
    old_cwd = os.getcwd()
    try:
        if xdmf_dir is not None:
            os.chdir(xdmf_dir)

        scene = GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()

        time_values = getattr(source, "TimestepValues", None)
        if not time_values:
            print("[INFO] Pas de série temporelle détectée (données stationnaires).")
            return

        if timestep_index == -1:
            target_time = time_values[-1]
        else:
            target_time = time_values[min(timestep_index, len(time_values) - 1)]

        scene.AnimationTime = target_time
        source.UpdatePipeline(target_time)
        print(f"[INFO] Positionné au temps t = {target_time}")
    finally:
        os.chdir(old_cwd)


def setup_view():
    """Configure la RenderView active : fond blanc, pas d'axes parasites."""
    view = GetActiveViewOrCreate("RenderView")
    view.Background = BACKGROUND_RGB
    view.UseColorPaletteForBackground = 0   # force le fond custom (pas le thème)
    view.OrientationAxesVisibility = 1 if SHOW_ORIENTATION_AXES else 0

    # Anti-aliasing pour des images propres à l'export
    try:
        view.EnableRayTracing = 0
    except Exception:
        pass

    return view


def apply_display_style(display):
    """Configure le style d'affichage de la représentation (edges, etc.)."""
    display.SetRepresentationType("Surface With Edges" if SHOW_EDGES else "Surface")
    if SHOW_CUBE_AXES:
        display.DataAxesGrid.GridAxesVisibility = 1
    else:
        display.DataAxesGrid.GridAxesVisibility = 0


def colorize_field(source, display, view, field_name, cfg):
    """Applique la colormap, les bornes et la légende propre pour un champ donné."""

    # Colorer par le champ (les champs vectoriels type U sont coloriés par leur magnitude)
    ColorBy(display, ("POINTS", field_name))

    # Récupère la fonction de transfert de couleur associée au champ
    ctf = GetColorTransferFunction(field_name)
    otf = GetOpacityTransferFunction(field_name)

    # Applique le preset choisi
    ctf.ApplyPreset(cfg["preset"], True)

    # Bornes : fixes si spécifiées, sinon rescale automatique sur les données
    if cfg["range"] is not None:
        ctf.RescaleTransferFunction(cfg["range"][0], cfg["range"][1])
        otf.RescaleTransferFunction(cfg["range"][0], cfg["range"][1])
    else:
        display.RescaleTransferFunctionToDataRange(False, True)

    # Echelle log si demandé
    ctf.UseLogScale = 1 if cfg["log"] else 0

    # Légende (scalar bar)
    display.SetScalarBarVisibility(view, True)
    sb = GetScalarBar(ctf, view)
    sb.Title = cfg.get("label", field_name)
    sb.ComponentTitle = ""
    sb.TitleColor = LEGEND_FONT_COLOR
    sb.LabelColor = LEGEND_FONT_COLOR
    sb.TitleFontSize = LEGEND_TITLE_FONT_SIZE
    sb.LabelFontSize = LEGEND_FONT_SIZE
    sb.Orientation = LEGEND_ORIENTATION
    sb.Position = LEGEND_POSITION
    sb.ScalarBarLength = 0.40 if LEGEND_ORIENTATION == "Horizontal" else 0.6
    sb.AutomaticLabelFormat = 0
    sb.LabelFormat = f"%-#6.{cfg['decimals']}f"
    sb.RangeLabelFormat = f"%-#6.{cfg['decimals']}f"
    sb.DrawBackground = 0
    sb.DrawScalarBarOutline = 0

    return ctf, sb


def export_current_view(view, outdir, field_name):
    """Exporte un screenshot PNG haute résolution du champ actuellement affiché."""
    os.makedirs(outdir, exist_ok=True)
    filepath = os.path.join(outdir, f"{field_name}.png")
    SaveScreenshot(
        filepath,
        view,
        ImageResolution=EXPORT_RESOLUTION,
        TransparentBackground=0,
    )
    print(f"[OK] Image exportée : {filepath}")


# =============================================================================
# PIPELINE PRINCIPAL
# =============================================================================

def run(input_path=None, outdir=DEFAULT_OUTDIR, fields=None, export=True, timestep=-1):
    fields_to_process = fields if fields else list(FIELDS.keys())

    source = get_active_source_or_load(input_path)

    xdmf_dir = None
    if input_path is not None and os.path.splitext(input_path)[1].lower() in (".xmf", ".xdmf"):
        xdmf_dir = os.path.dirname(os.path.abspath(input_path))

    goto_timestep(source, timestep, xdmf_dir=xdmf_dir)
    view = setup_view()
    display = Show(source, view)
    apply_display_style(display)

    Render()

    available_fields = [f.Name for f in source.PointData]
    print(f"[INFO] Champs disponibles dans le fichier : {available_fields}")

    for field_name in fields_to_process:
        if field_name not in FIELDS:
            print(f"[SKIP] '{field_name}' n'est pas défini dans FIELDS, on l'ignore.")
            continue
        if field_name not in available_fields:
            print(f"[SKIP] '{field_name}' absent des données chargées, on l'ignore.")
            continue

        cfg = FIELDS[field_name]
        print(f"[INFO] Traitement du champ : {field_name}")
        colorize_field(source, display, view, field_name, cfg)
        Render()

        if export:
            export_current_view(view, outdir, field_name)

    print("[DONE] Post-traitement terminé.")


# =============================================================================
# POINT D'ENTREE - utilisable en ligne de commande avec pvpython
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Post-traitement automatique ParaView - HYLON")
    parser.add_argument("--input", type=str, default=None,
                         help="Chemin vers le fichier de cas (.foam, .vtu, .pvd, ...). "
                              "Optionnel si une source est déjà chargée dans ParaView.")
    parser.add_argument("--outdir", type=str, default=DEFAULT_OUTDIR,
                         help="Dossier de sortie pour les images exportées.")
    parser.add_argument("--fields", type=str, nargs="*", default=None,
                         help="Liste de champs à traiter (sous-ensemble de FIELDS). "
                              "Par défaut : tous les champs définis dans FIELDS.")
    parser.add_argument("--no-export", action="store_true",
                         help="Ne pas exporter d'images, juste appliquer le style dans la vue.")
    parser.add_argument("--timestep", type=int, default=-1,
                         help="Index du pas de temps à afficher (-1 = dernier instant, défaut).")

    # pvpython passe parfois des arguments supplémentaires -> parse_known_args
    args, _ = parser.parse_known_args()

    run(
        input_path=args.input,
        outdir=args.outdir,
        fields=args.fields,
        export=not args.no_export,
        timestep=args.timestep,
    )