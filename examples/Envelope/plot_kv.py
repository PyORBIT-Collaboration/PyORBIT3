import os
import pathlib
import numpy as np
import matplotlib.pyplot as plt

from orbit.envelope import KVEnvelope

plt.style.use("style.mplstyle")


path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


envelope = KVEnvelope(
    eps_x=10.00e-06,
    eps_y=10.00e-06,
    mass=0.938,
    kin_energy=1.0,
    line_density=0.0,
    length=100.0,
)

x = envelope.sample(100_000)
xmax = np.std(x, axis=0) * 3.0
limits = list(zip(-xmax, xmax))

fig, axs = plt.subplots(ncols=4, nrows=4, figsize=(4, 4), constrained_layout=True)
for i in range(4):
    for j in range(4):
        ax = axs[i, j]
        ax.hist2d(x[:, j], x[:, i], bins=50, range=[limits[j], limits[i]], cmap="Blues")
for ax in axs.flat:
    ax.set_xticks([])
    ax.set_yticks([])
for i, dim in enumerate(["x", "xp", "y", "yp"]):
    axs[i, 0].set_ylabel(dim)
    axs[3, i].set_xlabel(dim)
fig.align_ylabels()

filename = "fig_corner.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()
