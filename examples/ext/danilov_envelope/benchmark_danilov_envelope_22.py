import argparse
import copy
import os
import pathlib
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

from orbit.danilov_envelope import DanilovEnvelope22


envelope = DanilovEnvelope22(
    intrinsic_emittance=1.0,
    eps_x_frac=0.5,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    intensity=0.0,
    mode=1,
    params=None
)

cov_matrix = envelope.cov()
print(cov_matrix)


x = envelope.sample(100_000, dist="gaussian")

xmax = 4.0 * np.std(x, axis=0)
limits = list(zip(-xmax, xmax))

fig, axs = plt.subplots(ncols=4, nrows=4, figsize=(6, 6), constrained_layout=True)
for i in range(4):
    for j in range(4):
        axis = (j, i)
        values, edges = np.histogramdd(x[:, axis], range=[limits[k] for k in axis], bins=50)
        axs[i, j].pcolormesh(edges[0], edges[1], values.T)
for ax in axs.flat:
    ax.set_xticks([])
    ax.set_yticks([])
plt.show()
