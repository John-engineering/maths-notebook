"""Shared plotting style for the notebook.

Palette: a validated categorical order (colour-blind safe on adjacent pairs).
Use SERIES[i] in fixed order; never cycle past 8 series.
"""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(os.path.dirname(HERE), "figures")

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
             "#256abf", "#184f95", "#0d366b"]
NEUTRAL = "#b9b8b3"


def apply():
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK_2,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "axes.axisbelow": True,
        "axes.prop_cycle": matplotlib.cycler(color=SERIES),
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "lines.linewidth": 2.0,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "font.size": 10,
        "figure.dpi": 110,
    })


def save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, bbox_inches="tight", dpi=130)
    plt.close(fig)
    return path


def plain_log(ax, axis="y"):
    """Log axis with plain-number tick labels (1, 10, 0.1 rather than 10^1)."""
    from matplotlib.ticker import FuncFormatter, NullFormatter
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    target = ax.yaxis if axis == "y" else ax.xaxis
    target.set_major_formatter(fmt)
    target.set_minor_formatter(NullFormatter())
