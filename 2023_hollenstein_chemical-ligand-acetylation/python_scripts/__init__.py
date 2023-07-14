from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def ticplot(
    x: np.ndarray,
    y: np.ndarray,
    ax: Optional[plt.Axes] = None,
    plotkwargs: Optional[dict] = None,
) -> plt.Axes:
    if ax is None:
        ax = plt.gca()
    if plotkwargs is None:
        plotkwargs = {"color": "#999999", "lw": 0.75}

    ax.plot(x, y, **plotkwargs)

    # Spines
    sns.despine(right=True, top=True, left=False, bottom=False)
    for spine in ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(0.5)
        spine.set_position(("outward", 5))

    # Adjust grind lines
    ax.yaxis.grid(
        which="both",
        linestyle="solid",
        # dashes=(4, 1.5),
        lw=0.5,
        alpha=1,
        color="#DDDDDD",
        zorder=0,
    )
    ax.xaxis.grid(False, which="both")

    # Tick and tick labels
    ax.tick_params(
        direction="out",
        length=1.5,
        width=0.5,
        colors="#000000",
        top=False,
        right=False,
        labelsize=6,
    )
    ax.yaxis.offsetText.set_fontsize(6)

    # Axis labels
    ax.set_xlabel("Retention Time (min)", fontsize=8)
    ax.set_ylabel("Total Ion Current", fontsize=8)
    return ax
