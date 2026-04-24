import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def calc_rms_ellipse_params(cov_matrix: np.ndarray) -> tuple[float, float, float]:
    """Return rms ellipse dimensions and orientation."""
    (i, j) = (0, 1)

    sii = cov_matrix[i, i]
    sjj = cov_matrix[j, j]
    sij = cov_matrix[i, j]

    angle = -0.5 * np.arctan2(2.0 * sij, sii - sjj)

    _sin = np.sin(angle)
    _cos = np.cos(angle)
    _sin2 = _sin**2
    _cos2 = _cos**2

    c1 = np.sqrt(abs(sii * _cos2 + sjj * _sin2 - 2 * sij * _sin * _cos))
    c2 = np.sqrt(abs(sii * _sin2 + sjj * _cos2 + 2 * sij * _sin * _cos))

    return (c1, c2, angle)


def plot_ellipse(
    r1: float = 1.0,
    r2: float = 1.0,
    angle: float = 0.0,
    center: tuple[float, float] = None,
    ax=None,
    **kws,
):
    kws.setdefault("fill", False)
    kws.setdefault("color", "black")
    kws.setdefault("lw", 1.25)

    center = (0.0, 0.0)

    d1 = r1 * 2.0
    d2 = r2 * 2.0
    angle = -np.degrees(angle)

    ax.add_patch(patches.Ellipse(center, d1, d2, angle=angle, **kws))
    return ax


def plot_rms_ellipse(
    cov_matrix: np.ndarray,
    level: float = 1.0,
    ax=None,
    **ellipse_kws,
):
    """Plot rms ellipse from 2 x 2 covariance matrix."""
    r1, r2, angle = calc_rms_ellipse_params(cov_matrix)
    plot_ellipse(r1 * level, r2 * level, angle=angle, ax=ax, **ellipse_kws)
    return ax


def plot_corner(
    particles: np.ndarray,
    limits: list[tuple[float, float]] = None,
    bins: int = 64,
    labels: list[str] = None,
    blur: float = None,
) -> tuple:
    """Generate corner plot."""
    ndim = particles.shape[1]

    if limits is None:
        xmax = np.max(particles, axis=0)
        xmin = np.min(particles, axis=0)
        limits = list(zip(xmin, xmax))

    if labels is None:
        labels = ndim * [""]

    fig, axs = plt.subplots(
        ncols=ndim, nrows=ndim, sharex=None, sharey=None, figsize=(8, 8)
    )
    for i in range(ndim):
        for j in range(ndim):
            axis = (j, i)
            ax = axs[i, j]
            if i > j:
                values, edges = np.histogramdd(
                    particles[:, axis], bins=bins, range=[limits[k] for k in axis]
                )
                if blur:
                    values = scipy.ndimage.gaussian_filter(values, sigma=blur)
                ax.pcolormesh(
                    edges[0],
                    edges[1],
                    values.T,
                    linewidth=0.0,
                    rasterized=True,
                    shading="auto",
                )
            elif i == j:
                values, edges = np.histogram(
                    particles[:, i], bins=bins, range=limits[i]
                )
                if blur:
                    values = scipy.ndimage.gaussian_filter(values, sigma=blur)
                ax.stairs(values, edges, lw=1.5, color="black")
            else:
                ax.axis("off")

    for i in range(0, ndim - 1):
        for j in range(0, ndim):
            axs[i, j].set_xticklabels([])
    for i in range(0, ndim):
        for j in range(1, ndim):
            axs[i, j].set_yticklabels([])

    for ax in axs.flat:
        for loc in ["top", "right"]:
            ax.spines[loc].set_visible(False)

    for i, label in enumerate(labels):
        axs[-1, i].set_xlabel(label)
    for i, label in enumerate(labels[1:], start=1):
        axs[i, 0].set_ylabel(label)

    axs[0, 0].set_yticklabels([])
    axs[0, 0].set_ylabel(None)

    fig.align_ylabels()
    fig.align_xlabels()

    return fig, axs
