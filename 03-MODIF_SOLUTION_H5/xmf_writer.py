"""
xmf_writer.py
=============
Génération du fichier XMF (XDMF2) pour visualisation Paraview.
À intégrer dans modify_cfd_solution.py ou utiliser en standalone.

Deux modes :
  - write_xmf_dual()   : XMF pointant vers mesh.h5 + sol.h5 séparément
  - write_xmf_single() : XMF pointant vers un sol.h5 qui contient déjà les coords
"""

import h5py
import numpy as np
import os


# ═══════════════════════════════════════════════════════════════════
#  DÉTECTION AUTOMATIQUE DE LA CONNECTIVITÉ DANS LE MESH
# ═══════════════════════════════════════════════════════════════════

# Types d'éléments XDMF2 et leur topologie
# https://www.xdmf.org/index.php/XDMF_Model_and_Format
XDMF_TOPOLOGY = {
    "tetra"   : ("Tetrahedron",  4),
    "tetra10" : ("Tet_10",      10),
    "hexa"    : ("Hexahedron",   8),
    "tri"     : ("Triangle",     3),
    "quad"    : ("Quadrilateral",4),
    "wedge"   : ("Wedge",        6),
    "pyramid" : ("Pyramid",      5),
    # YALES2 maillages mixtes
    "mixed"   : ("Mixed",        None),
}

# Chemins candidats pour la connectivité dans les mesh
CONNECTIVITY_CANDIDATES = [
    # (groupe, dataset, type_element)
    ("Connectivity", "tet->node", "tetra"),
    ("Topology",          "Tetrahedra",     "tetra"),
    ("Topology",          "Hexahedra",      "hexa"),
    ("Topology",          "connectivity",   "mixed"),
    ("",                  "Tetrahedra",     "tetra"),
    ("",                  "connectivity",   "mixed"),
    ("Mesh/Topology",     "Tetrahedra",     "tetra"),
    ("Mesh",              "Tetrahedra",     "tetra"),
    ("Mesh",              "connectivity",   "mixed"),
    ("Elements",          "Tetrahedra",     "tetra"),
    ("Elements/Tetrahedra","connectivity",  "tetra"),
]

# Chemins candidats pour les coordonnées
COORD_CANDIDATES = [
    ("Coordinates", "x", "y", "z"),
    ("Coordinates", "X", "Y", "Z"),
    ("Nodes",       "x", "y", "z"),
    ("Nodes",       "X", "Y", "Z"),
    ("Mesh/Coordinates", "x", "y", "z"),
    ("",            "x", "y", "z"),
    ("",            "X", "Y", "Z"),
]


def probe_mesh_structure(mesh_file: str) -> dict:
    """
    Inspecte le fichier mesh et retourne un dict avec :
      - coord_group, coord_x/y/z  : chemins des coordonnées
      - conn_group, conn_dataset   : chemin de la connectivité
      - elem_type                  : type d'élément XDMF
      - n_nodes, n_elems, n_pts_per_elem
    """
    info = {}
    with h5py.File(mesh_file, "r") as f:

        # ── Coordonnées ──────────────────────────────────────────
        for (grp, dx, dy, dz) in COORD_CANDIDATES:
            try:
                base = f[grp] if grp else f
                x = base[dx]
                info["coord_group"] = grp
                info["coord_x"]     = dx
                info["coord_y"]     = dy
                info["coord_z"]     = dz
                info["n_nodes"]     = len(x)
                info["coord_dtype_prec"] = 8 if x.dtype == np.float64 else 4
                print(f"  ✓ Coords dans '{grp or '/'}' : {dx},{dy},{dz}  "
                      f"({info['n_nodes']:,} nœuds)")
                break
            except (KeyError, TypeError):
                continue

        if "n_nodes" not in info:
            print("  ⚠ Coordonnées non trouvées automatiquement.")
            _dump_structure(f)
            raise KeyError("Coordonnées introuvables. Voir la structure ci-dessus.")

        # ── Connectivité ─────────────────────────────────────────
        for (grp, ds, etype) in CONNECTIVITY_CANDIDATES:
            try:
                base = f[grp] if grp else f
                conn = base[ds]
                shape = conn.shape
                info["conn_group"]      = grp
                info["conn_dataset"]    = ds
                info["elem_type"]       = etype
                if etype != "mixed":
                    n_pts = XDMF_TOPOLOGY[etype][1]

                    if len(shape) == 1:
                        # connectivité stockée à plat
                        info["n_elems"] = shape[0] // n_pts
                    else:
                        info["n_elems"] = shape[0]

                    info["n_pts_per_elem"] = n_pts
                    conn = base[ds]
                    shape = conn.shape

                    info["conn_shape"] = shape
                    info["conn_group"] = grp
                    info["conn_dataset"] = ds
                    info["elem_type"] = etype
                else:
                    info["n_elems"]          = shape[0]   # à affiner si mixed
                    info["n_pts_per_elem"]   = None
                info["conn_dtype_prec"] = 8 if conn.dtype in (np.int64,) else 4
                print(f"  ✓ Connectivité dans '{grp or '/'}' / '{ds}' : "
                      f"shape={shape}  type={etype}")
                break
            except (KeyError, TypeError):
                continue

        if "conn_group" not in info:
            print("  ⚠ Connectivité non trouvée automatiquement.")
            _dump_structure(f)
            raise KeyError("Connectivité introuvable. Voir la structure ci-dessus.")

    return info


def _dump_structure(f):
    """Affiche la structure HDF5 pour diagnostic."""
    print("\n  ── Structure du fichier ──")
    def show(name, obj):
        tag = "DS" if isinstance(obj, h5py.Dataset) else "GR"
        shp = f" {obj.shape}" if isinstance(obj, h5py.Dataset) else ""
        print(f"    [{tag}] {name}{shp}")
    f.visititems(show)


# ═══════════════════════════════════════════════════════════════════
#  STRATÉGIE A : XMF pointant vers mesh.h5 + sol.h5 séparément
# ═══════════════════════════════════════════════════════════════════

def write_xmf_dual(
    xmf_file: str,
    mesh_file: str,
    sol_file:  str,
    mesh_info: dict,
    time_value: float = 0.0,
):
    """
    Génère un XMF qui référence :
      - les coordonnées et la connectivité depuis mesh_file
      - tous les champs scalaires depuis sol_file

    Les chemins dans le XMF sont RELATIFS au répertoire du XMF.
    Mets les 3 fichiers dans le même dossier.
    """
    # Chemins relatifs
    xmf_dir  = os.path.dirname(os.path.abspath(xmf_file))
    mesh_rel = os.path.relpath(mesh_file, xmf_dir)
    sol_rel  = os.path.relpath(sol_file,  xmf_dir)

    nnode = mesh_info["n_nodes"]
    nelems = mesh_info["n_elems"]
    etype, npts = XDMF_TOPOLOGY[mesh_info["elem_type"]]

    # Chemin HDF5 de la connectivité
    conn_h5_path = (
        f"{mesh_info['conn_group']}/{mesh_info['conn_dataset']}"
        if mesh_info["conn_group"]
        else mesh_info["conn_dataset"]
    )
    conn_prec = mesh_info["conn_dtype_prec"]

    # Chemin HDF5 des coordonnées
    cg = mesh_info["coord_group"]
    cx = mesh_info["coord_x"]
    cy = mesh_info["coord_y"]
    cz = mesh_info["coord_z"]
    coord_prec = mesh_info["coord_dtype_prec"]

    # Collecter tous les champs depuis sol_file
    fields = _collect_fields(sol_file, nnode)

    lines = []
    lines.append('<?xml version="1.0" ?>')
    lines.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    lines.append('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">')
    lines.append('  <Domain>')
    lines.append('    <Grid Name="mesh" GridType="Uniform">')
    lines.append(f'      <Time Value="{time_value:.8e}"/>')

    # ── Topologie ────────────────────────────────────────────────
    raw_conn_size = mesh_info["conn_shape"][0]

    lines.append(f'      <Topology Type="Tetrahedron" NumberOfElements="{nelems}">')
    lines.append(f'        <DataItem ItemType="Function" ' f'Dimensions="{raw_conn_size}" ' f'Function="$0 - 1">'    )
    lines.append(f'          <DataItem Format="HDF" ' f'DataType="Int" ' f'Dimensions="{raw_conn_size}">' )
    lines.append(f'           {mesh_rel}:/{conn_h5_path}')
    lines.append('          </DataItem>')
    lines.append('        </DataItem>')
    lines.append('      </Topology>')

    # ── Géométrie (coordonnées entrelacées X Y Z) ─────────────────
    # XDMF supporte XYZ séparés (X_Y_Z) ou entrelacés (XYZ)
    # On utilise X_Y_Z pour pointer directement vers les datasets séparés
    lines.append('      <Geometry Type="X_Y_Z">')
    for axis, ds in [(cx, cx), (cy, cy), (cz, cz)]:
        h5path = f"{cg}/{ds}" if cg else ds
        lines.append(f'        <DataItem Dimensions="{nnode}" '
                     f'NumberType="Float" Precision="{coord_prec}" Format="HDF">')
        lines.append(f'          {mesh_rel}:/{h5path}')
        lines.append( '        </DataItem>')
    lines.append( '      </Geometry>')

    # ── Champs ───────────────────────────────────────────────────
    for (attr_name, h5path, prec) in fields:
        lines.append(f'      <Attribute Name="{attr_name}" 'f'Center="Node" AttributeType="Scalar">')
        lines.append(f'        <DataItem Dimensions="{nnode}" '
                     f'NumberType="Float" Precision="{prec}" Format="HDF">')
        lines.append(f'          {sol_rel}:/{h5path}')
        lines.append( '        </DataItem>')
        lines.append( '      </Attribute>')

    lines.append('    </Grid>')
    lines.append('  </Domain>')
    lines.append('</Xdmf>')

    with open(xmf_file, "w") as fout:
        fout.write("\n".join(lines) + "\n")

    print(f"  ✓ XMF écrit : {xmf_file}")
    print(f"    → {len(fields)} champs référencés")


# ═══════════════════════════════════════════════════════════════════
#  STRATÉGIE B : coords copiées dans sol.h5 → XMF fichier unique
# ═══════════════════════════════════════════════════════════════════

def embed_coords_in_solution(
    mesh_file: str,
    sol_file:  str,
    mesh_info: dict,
):
    """
    Copie les coordonnées et la connectivité du mesh dans le fichier solution.
    Après cela, on n'a besoin que du sol.h5 pour visualiser.
    ⚠ Augmente la taille du fichier sol d'environ 75 Mo pour 3.2M nœuds.
    """
    cg = mesh_info["coord_group"]
    cx = mesh_info["coord_x"]
    cy = mesh_info["coord_y"]
    cz = mesh_info["coord_z"]
    conn_grp = mesh_info["conn_group"]
    conn_ds  = mesh_info["conn_dataset"]

    with h5py.File(mesh_file, "r") as fm, h5py.File(sol_file, "r+") as fs:
        # Coords
        base_m = fm[cg] if cg else fm
        for ds_name in [cx, cy, cz]:
            dest = f"Mesh/Coordinates/{ds_name}"
            if dest not in fs:
                fs.create_dataset(dest, data=base_m[ds_name][:])
        # Connectivité
        base_c = fm[conn_grp] if conn_grp else fm
        dest_c = f"Mesh/Topology/{conn_ds}"
        if dest_c not in fs:
            fs.create_dataset(dest_c, data=base_c[conn_ds][:])

    print(f"  ✓ Coords + connectivité copiées dans {sol_file}")


def write_xmf_single(
    xmf_file:  str,
    sol_file:  str,
    mesh_info: dict,
    time_value: float = 0.0,
):
    """
    Génère un XMF autonome : toutes les données dans sol_file.
    (Nécessite d'avoir appelé embed_coords_in_solution() avant.)
    """
    xmf_dir = os.path.dirname(os.path.abspath(xmf_file))
    sol_rel  = os.path.relpath(sol_file, xmf_dir)

    cx = mesh_info["coord_x"]
    cy = mesh_info["coord_y"]
    cz = mesh_info["coord_z"]
    conn_ds  = mesh_info["conn_dataset"]
    etype, npts = XDMF_TOPOLOGY[mesh_info["elem_type"]]
    nnode  = mesh_info["n_nodes"]
    nelems = mesh_info["n_elems"]
    conn_prec  = mesh_info["conn_dtype_prec"]
    coord_prec = mesh_info["coord_dtype_prec"]

    fields = _collect_fields(sol_file, nnode)

    lines = []
    lines.append('<?xml version="1.0" ?>')
    lines.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    lines.append('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">')
    lines.append('  <Domain>')
    lines.append('    <Grid Collection="Tetrahedron_Mesh" ''Name="solution-Tetrahedron">')
    lines.append(f'      <Time Value="{time_value:.8e}"/>')

    lines.append(f'      <Topology TopologyType="{etype}" NumberOfElements="{nelems}">')
    lines.append(f'        <DataItem Dimensions="{nelems} {npts}" '
                 f'NumberType="Int" Precision="{conn_prec}" Format="HDF">')
    lines.append(f'          {sol_rel}:/Mesh/Topology/{conn_ds}')
    lines.append( '        </DataItem>')
    lines.append( '      </Topology>')

    lines.append( '      <Geometry GeometryType="X_Y_Z">')
    for axis in [cx, cy, cz]:
        lines.append(f'        <DataItem Dimensions="{nnode}" '
                     f'NumberType="Float" Precision="{coord_prec}" Format="HDF">')
        lines.append(f'          {sol_rel}:/Mesh/Coordinates/{axis}')
        lines.append( '        </DataItem>')
    lines.append( '      </Geometry>')

    for (attr_name, h5path, prec) in fields:
        lines.append(f'      <Attribute Name="{attr_name}" '
                     f'AttributeType="Scalar" Center="Node">')
        lines.append(f'        <DataItem Dimensions="{nnode}" '
                     f'NumberType="Float" Precision="{prec}" Format="HDF">')
        lines.append(f'          {sol_rel}:/{h5path}')
        lines.append( '        </DataItem>')
        lines.append( '      </Attribute>')

    lines.append('    </Grid>')
    lines.append('  </Domain>')
    lines.append('</Xdmf>')

    with open(xmf_file, "w") as fout:
        fout.write("\n".join(lines) + "\n")

    print(f"  ✓ XMF autonome écrit : {xmf_file}")


# ═══════════════════════════════════════════════════════════════════
#  COLLECTE DES CHAMPS DEPUIS LE FICHIER SOLUTION
# ═══════════════════════════════════════════════════════════════════

def _collect_fields(sol_file: str, nnode: int) -> list:
    """
    Retourne la liste de (nom_affichage, chemin_hdf5, precision)
    pour tous les champs scalaires nodaux du fichier solution.
    Calcule aussi Y_k = rhoY_k / rho et les ajoute comme champs dérivés.
    """
    fields = []

    with h5py.File(sol_file, "r") as f:
        rho_path = "GaseousPhase/rho"

        # ── Variables conservatives brutes ───────────────────────
        conservative_skip = {"rho"}  # rho ajouté séparément
        for ds_name in f["GaseousPhase"]:
            path = f"GaseousPhase/{ds_name}"
            d = f[path]
            if d.shape == (nnode,):
                prec = 8 if d.dtype == np.float64 else 4
                fields.append((ds_name, path, prec))

        # ── Fractions massiques Y_k = rhoY_k / rho ───────────────
        # On ne peut pas les stocker directement dans le XMF (calcul)
        # → on les PRÉ-CALCULE et les écrit dans le sol.h5 comme dataset temporaire
        species = list(f["RhoSpecies"].keys())
        rho = f[rho_path][:]

        # Vérifier si les Y_k sont déjà présents (run précédent)
        for sp in species:
            y_path = f"MassFractions/Y_{sp}"
            if y_path in f:
                d = f[y_path]
                prec = 8 if d.dtype == np.float64 else 4
                fields.append((f"Y_{sp}", y_path, prec))
            else:
                # Sera écrit juste après (voir write_mass_fractions)
                fields.append((f"Y_{sp}", y_path, 4))

        # ── Additionals ──────────────────────────────────────────
        for ds_name in f["Additionals"]:
            path = f"Additionals/{ds_name}"
            d = f[path]
            if d.shape == (nnode,):
                prec = 8 if d.dtype == np.float64 else 4
                fields.append((ds_name, path, prec))

    return fields


def write_mass_fractions(sol_file: str):
    """
    Calcule et écrit Y_k = rhoY_k / rho dans le groupe MassFractions/
    du fichier solution, pour que le XMF puisse y pointer directement.
    """
    with h5py.File(sol_file, "r+") as f:
        rho     = f["GaseousPhase/rho"][:]
        species = list(f["RhoSpecies"].keys())

        if "MassFractions" not in f:
            f.create_group("MassFractions")

        for sp in species:
            y_path = f"MassFractions/Y_{sp}"
            rhoY   = f[f"RhoSpecies/{sp}"][:]
            Y      = (rhoY / rho).astype(np.float32)

            if y_path in f:
                f[y_path][:] = Y           # mise à jour
            else:
                f.create_dataset(y_path, data=Y, dtype=np.float32)

        print(f"  ✓ Fractions massiques Y_k écrites dans MassFractions/ "
              f"({len(species)} espèces)")


# ═══════════════════════════════════════════════════════════════════
#  POINT D'ENTRÉE PRINCIPAL (à appeler depuis modify_cfd_solution.py)
# ═══════════════════════════════════════════════════════════════════

def generate_xmf(
    sol_file:   str,
    mesh_file:  str = None,
    xmf_file:   str = None,
    time_value: float = 0.0,
    mode:       str = "dual",   # "dual" ou "single"
):
    """
    Fonction principale à appeler après la modification de la solution.

    Parameters
    ----------
    sol_file   : fichier solution modifié (.h5)
    mesh_file  : fichier maillage (.h5)  — requis si mode="dual"
    xmf_file   : chemin de sortie du XMF  [défaut : même base que sol_file]
    time_value : temps physique à encoder dans le XMF
    mode       : "dual"   → XMF pointe vers mesh.h5 + sol.h5
                 "single" → coordonnées copiées dans sol.h5 (autonome)
    """
    if xmf_file is None:
        base     = os.path.splitext(sol_file)[0]
        xmf_file = base + ".xmf"

    print(f"\n── Génération du fichier XMF ({mode}) ──")

    # Écrire les fractions massiques Y_k dans le sol.h5
    write_mass_fractions(sol_file)

    if mode == "dual":
        if mesh_file is None:
            raise ValueError("mesh_file requis pour le mode 'dual'.")
        print("  Inspection du maillage …")
        mesh_info = probe_mesh_structure(mesh_file)
        write_xmf_dual(xmf_file, mesh_file, sol_file, mesh_info, time_value)

    elif mode == "single":
        if mesh_file is None:
            raise ValueError("mesh_file requis pour copier les coordonnées.")
        print("  Inspection du maillage …")
        mesh_info = probe_mesh_structure(mesh_file)
        embed_coords_in_solution(mesh_file, sol_file, mesh_info)
        write_xmf_single(xmf_file, sol_file, mesh_info, time_value)

    else:
        raise ValueError(f"mode inconnu : '{mode}'. Choisir 'dual' ou 'single'.")

    return xmf_file


# ── Test standalone ──────────────────────────────────────────────
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage : python xmf_writer.py sol.h5 mesh.h5 [dual|single]")
        sys.exit(1)
    sol  = sys.argv[1]
    mesh = sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "dual"
    generate_xmf(sol, mesh, mode=mode)