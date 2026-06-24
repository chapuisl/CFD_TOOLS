import Graphic_Configuration as GC

"""
# ===================================================================================================================
#  Import library
# ===================================================================================================================
"""
from matplotlib.ticker import LogLocator, LogFormatterSciNotation, NullFormatter, FuncFormatter
from matplotlib.legend_handler import HandlerTuple
import matplotlib.pyplot as plt
import numpy as np
import os

"""
# ===================================================================================================================
#  Function
# ===================================================================================================================
"""

def _to_list(x):
    """Wrap a single array/list into a list-of-lists."""
    if isinstance(x, (list, tuple)) and len(x) > 0 and isinstance(x[0], (list, tuple, np.ndarray)):
        return x
    return [x]


def _nice_ticks(vmin, vmax, n_target=5):
    span = vmax - vmin
    if span == 0:
        return np.array([vmin]), vmin, vmax

    best = None
    best_score = np.inf

    mag = 10 ** np.floor(np.log10(span))
    for exp_offset in [-1, 0, 1]:
        m = mag * (10 ** exp_offset)
        for factor in [1, 2, 2.5, 5]:
            step = factor * m
            new_vmin = np.floor(vmin / step) * step
            new_vmax = np.ceil(vmax / step) * step
            ticks = np.arange(new_vmin, new_vmax + step * 0.5, step)
            n = len(ticks)
            overshoot = ((new_vmax - vmax) + (vmin - new_vmin)) / span
            score = abs(n - n_target) * 2 + overshoot
            if score < best_score:
                best_score = score
                best = (ticks, new_vmin, new_vmax)

    return best


def _set_axis_ticks(ax, scale_x, scale_y, xmin, xmax, ymin, ymax):
    for axis, scale, vmin, vmax in [('x', scale_x, xmin, xmax),
                                     ('y', scale_y, ymin, ymax)]:
        if scale == 'log':
            continue
        setter    = ax.set_xlim   if axis == 'x' else ax.set_ylim
        set_ticks = ax.set_xticks if axis == 'x' else ax.set_yticks
        fmt_axis  = ax.xaxis      if axis == 'x' else ax.yaxis

        ticks, new_min, new_max = _nice_ticks(vmin, vmax, n_target=5)
        setter(new_min, new_max)
        set_ticks(ticks)

        if abs(new_max) > 1e4 or abs(new_max) < 1e-2 or (new_max != 0 and new_min != 0 and abs(new_min / new_max) < 1e-2):
            fmt_axis.set_major_formatter(FuncFormatter(
                lambda x, _: f'{x:.2e}' if x != 0 else '0'
            ))
        else:
            step = ticks[1] - ticks[0] if len(ticks) > 1 else 1
            dec = max(0, -int(np.floor(np.log10(abs(step)))) if step != 0 else 0)
            fmt_axis.set_major_formatter(FuncFormatter(
                lambda x, _, d=dec: f'{x:.{d}f}'
            ))


def _set_log_axis(ax, axis):
    a = ax.xaxis if axis == 'x' else ax.yaxis
    a.set_major_locator(LogLocator(base=10, numticks=10))
    a.set_major_formatter(LogFormatterSciNotation(base=10))
    a.set_minor_locator(LogLocator(base=10, subs='auto', numticks=10))
    a.set_minor_formatter(LogFormatterSciNotation(base=10, minor_thresholds=(2, 0.5)))
    tick_axis = 'x' if axis == 'x' else 'y'
    ax.tick_params(axis=tick_axis, which='minor', labelsize=14, rotation=45)


def _add_ref_lines(ax, line_value, line_orientation, type_x_scale, type_y_scale):
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_range, y_range = x_max - x_min, y_max - y_min

    offset_h = offset_v = 0
    for val, ori in zip(line_value, line_orientation):
        if ori == 'H':
            ax.axhline(y=val, color='red', linestyle='--', linewidth=2)
            x_text = x_max + 0.2 * x_range
            y_text = val + offset_h * 0.1 * y_range
            ax.annotate(f'{val:.2e}', xy=(x_max, val), xytext=(x_text, y_text),
                        color='red', va='center', ha='left',
                        arrowprops=dict(arrowstyle='-', color='red', lw=0.8),
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
            offset_h += 1
        elif ori == 'V':
            ax.axvline(x=val, color='b', linestyle='--', linewidth=2)
            x_text = val + offset_v * 0.1 * x_range
            y_text = y_max + 0.2 * y_range
            ax.annotate(f'{val:.2e}', xy=(val, y_max), xytext=(x_text, y_text),
                        color='b', va='bottom', ha='right', rotation=90,
                        arrowprops=dict(arrowstyle='-', color='b', lw=0.8),
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
            offset_v += 1


def _auto_rotate_xlabels(ax, fig, threshold=0.9):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    formatter = ax.xaxis.get_major_formatter()
    ticks = ax.get_xticks()
    tick_strings = [formatter(t, i) for i, t in enumerate(ticks)]
    tick_strings = [s for s in tick_strings if s.strip()]
    if not tick_strings:
        return
    tmp_text = ax.text(0, 0, max(tick_strings, key=len),
                       fontsize=plt.rcParams['xtick.labelsize'])
    bbox = tmp_text.get_window_extent(renderer=renderer)
    max_label_w_px = bbox.width
    tmp_text.remove()
    n_labels = len(tick_strings)
    axis_w_px = ax.get_window_extent(renderer=renderer).width
    space_per_tick_px = axis_w_px / n_labels
    angle = 45 if max_label_w_px > threshold * space_per_tick_px else 0
    ha    = 'right' if angle == 45 else 'center'
    ax.tick_params(axis='x', labelrotation=angle)
    plt.setp(ax.get_xticklabels(), ha=ha)


# ───────────────────────────────────────────────────────────────
#  Fonction principale
# ───────────────────────────────────────────────────────────────

def plot_evolution(
    t, data,
    labels=None, colors=None, styles=None,
    xlabel=None, ylabel=None, figsize=None, title="", legend_loc="best",
    x_limit_left=None, x_limit_right=None, y_limit_bot=None, y_limit_top=None,
    type_x_scale='linear', type_y_scale='linear',
    line_value=None, line_orientation=None,
    secondary_data=None, secondary_label=None, secondary_color=None, secondary_ylabel=None,
    secondary_axis_offset=70,
    marker=None, plot_fig=False, save_fig=False, save_path=None, name_fig=None,):

    def _is_scalar_series(x):
        """True si x est une série 1-D (une seule courbe)."""
        if isinstance(x, np.ndarray):
            return x.ndim == 1
        if isinstance(x, (list, tuple)):
            if len(x) == 0:
                return True
            return not isinstance(x[0], (list, tuple, np.ndarray))
        return True   # scalaire isolé → série dégénérée

    def _normalize_to_axis_list(sd):
        """
        Retourne une liste de longueur N_axes.
        Chaque élément est lui-même une liste de séries (pour _to_list).
        """
        if _is_scalar_series(sd):
            # Cas (A) : une seule courbe, un seul axe
            return [sd]

        # sd est une liste dont les éléments sont soit des scalaires soit des arrays
        # On regarde si TOUS les éléments sont des séries scalaires
        all_scalar = all(_is_scalar_series(s) for s in sd)

        if all_scalar:
            # Cas (B) : plusieurs courbes sur un seul axe
            return [sd]

        # Cas (C) : chaque élément = un axe (qui peut lui-même être (A) ou (B))
        # On ré-applique récursivement pour normaliser chaque axe en liste de séries
        result = []
        for elem in sd:
            if _is_scalar_series(elem):
                result.append(elem)           # axe avec 1 courbe
            else:
                result.append(elem)           # axe avec N courbes (_to_list s'en chargera)
        return result

    if secondary_data is not None:
        sec_data_list = _normalize_to_axis_list(secondary_data)
        N = len(sec_data_list)

        # ── secondary_label ──
        if secondary_label is None:
            sec_label_list = [None] * N
        elif isinstance(secondary_label, str):
            sec_label_list = [secondary_label] + [None] * (N - 1)
        else:
            # liste de labels (un par axe) ou liste de listes (plusieurs courbes par axe)
            sec_label_list = list(secondary_label) + [None] * max(0, N - len(secondary_label))

        # ── secondary_color ──
        if secondary_color is None:
            sec_color_list = [None] * N
        elif not isinstance(secondary_color, (list, tuple)):
            sec_color_list = [secondary_color] + [None] * (N - 1)
        else:
            sec_color_list = list(secondary_color) + [None] * max(0, N - len(secondary_color))

        # ── secondary_ylabel ──
        if secondary_ylabel is None:
            sec_ylabel_list = [None] * N
        elif isinstance(secondary_ylabel, str):
            sec_ylabel_list = [secondary_ylabel] + [None] * (N - 1)
        else:
            sec_ylabel_list = list(secondary_ylabel) + [None] * max(0, N - len(secondary_ylabel))
    else:
        sec_data_list = []
        sec_label_list = sec_color_list = sec_ylabel_list = []

    # Palette de couleurs par défaut pour les axes secondaires sans couleur définie
    _default_sec_colors = ['darkred', 'darkorange', 'darkgreen', 'darkcyan', 'darkmagenta']

    # ── Création de la figure ──────────────────────────────────────────────────
    # Élargit automatiquement la figure si plusieurs axes secondaires
    if figsize is None:
        base_w = 10
        extra_w = max(0, len(sec_data_list) - 1) * 1.2   # ~1.2 pouce par axe sup.
        figsize = (base_w + extra_w, 6)

    fig, ax1 = plt.subplots(figsize=figsize)
    ax1.grid(True, which="both")
    ax1.set_xscale(type_x_scale)
    ax1.set_yscale(type_y_scale)

    # ── Tracé des courbes principales ─────────────────────────────────────────
    data_list = _to_list(data)
    t_list    = _to_list(t)
    if len(t_list) == 1:
        t_list = t_list * len(data_list)

    primary_lines = []
    for i, (ti, di) in enumerate(zip(t_list, data_list)):
        l, = ax1.plot(ti, di,
                 color=colors[i] if colors and i < len(colors) else '0.5',
                 linestyle=styles[i] if styles and i < len(styles) else '-',
                 label=labels[i] if labels and i < len(labels) else None,
                 marker=marker)
        primary_lines.append(l)

    if type_x_scale == 'log':
        _set_log_axis(ax1, 'x')
    if type_y_scale == 'log':
        _set_log_axis(ax1, 'y')

    # ── Bornes réelles des données ────────────────────────────────────────────
    all_x = np.concatenate([np.asarray(ti).ravel() for ti in t_list])
    all_y = np.concatenate([np.asarray(di).ravel() for di in data_list])
    raw_xmin = x_limit_left  if x_limit_left  is not None else float(np.nanmin(all_x))
    raw_xmax = x_limit_right if x_limit_right is not None else float(np.nanmax(all_x))
    raw_ymin = y_limit_bot   if y_limit_bot   is not None else float(np.nanmin(all_y))
    raw_ymax = y_limit_top   if y_limit_top   is not None else float(np.nanmax(all_y))

    _set_axis_ticks(ax1, type_x_scale, type_y_scale, raw_xmin, raw_xmax, raw_ymin, raw_ymax)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    _auto_rotate_xlabels(ax1, fig)

    # ── Lignes de référence ───────────────────────────────────────────────────
    if line_value is not None and line_orientation is not None:
        _add_ref_lines(ax1, line_value, line_orientation, type_x_scale, type_y_scale)

    # ── Axes secondaires (0, 1 ou N) ─────────────────────────────────────────
    secondary_axes = []
    secondary_lines = []

    for k, (sec_data, sec_label, sec_color, sec_ylabel) in enumerate(
            zip(sec_data_list, sec_label_list, sec_color_list, sec_ylabel_list)):

        # Couleur par défaut
        if sec_color is None:
            if k == 0:
                sec_color = GC.Purple_Black()
            else:
                sec_color = _default_sec_colors[k % len(_default_sec_colors)]

        ax_sec = ax1.twinx()

        # Décalage de la spine droite pour les axes k >= 1
        if k > 0:
            ax_sec.spines["right"].set_position(("outward", secondary_axis_offset * k))

        # Tracé de la série secondaire
        sec_series_list = _to_list(sec_data)
        t_sec = t_list if len(t_list) == len(sec_series_list) else [t_list[0]] * len(sec_series_list)

        for j, (ti, di) in enumerate(zip(t_sec, sec_series_list)):
            lbl = sec_label if j == 0 else None   # un seul label par axe
            l, = ax_sec.plot(ti, di,
                             color=sec_color,
                             linestyle='--',
                             label=lbl,
                             marker=marker)
            secondary_lines.append(l)

        # Titre et couleur de l'axe Y secondaire
        ylabel_sec = sec_ylabel if sec_ylabel is not None else (sec_label or "")
        ax_sec.set_ylabel(ylabel_sec, color=sec_color)
        ax_sec.tick_params(axis='y', labelcolor=sec_color)

        # Ticks propres pour cet axe
        all_sy = np.concatenate([np.asarray(di).ravel() for di in sec_series_list])
        all_sx = np.concatenate([np.asarray(ti).ravel() for ti in t_sec])
        _set_axis_ticks(ax_sec, 'linear', 'linear',
                        float(np.nanmin(all_sx)), float(np.nanmax(all_sx)),
                        float(np.nanmin(all_sy)), float(np.nanmax(all_sy)))

        secondary_axes.append(ax_sec)

    # ── Légende combinée sur ax1 ──────────────────────────────────────────────
    all_lines  = primary_lines + secondary_lines
    all_labels = [l.get_label() for l in all_lines]
    # Filtre les labels vides ou auto-générés (_line…)
    filtered = [(l, lb) for l, lb in zip(all_lines, all_labels)
                if lb and not lb.startswith('_')]
    if filtered:
        lns, lbs = zip(*filtered)
        ax1.legend(lns, lbs, loc=legend_loc)

    # ── Ajustement de la marge droite si plusieurs axes secondaires ──────────
    if len(secondary_axes) > 1:
        # right doit laisser de la place pour les axes décalés
        right_margin = 1 - 0.06 * len(secondary_axes)
        fig.subplots_adjust(right=max(0.60, right_margin))

    # ── Sauvegarde ────────────────────────────────────────────────────────────
    if save_fig and save_path:
        os.makedirs(save_path, exist_ok=True)
        fname = name_fig.replace(" ", "_").replace("/", "_").replace("\\", "_")
        plt.savefig(os.path.join(save_path, f"{fname}.png"), dpi=300, bbox_inches='tight')
        plt.close()

    if plot_fig:
        plt.show()