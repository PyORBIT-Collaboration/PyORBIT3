import numpy as np
import matplotlib.pyplot as plt

from orbit.danilov_envelope import DanilovEnvelope22

plt.style.use("style.mplstyle")


envelope = DanilovEnvelope22(
    intrinsic_emittance=20.00e-06,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    intensity=0.0,
    mode=1,
    params=None,
)

x = envelope.sample(100_000)
xmax = np.std(x, axis=0) * 3.0
limits = list(zip(-xmax, xmax))

fig, axs = plt.subplots(ncols=4, nrows=4, figsize=(4, 4), constrained_layout=True)
for i in range(4):
    for j in range(4):
        ax = axs[i, j]
        ax.hist2d(x[:, j], x[:, i], bins=50, range=[limits[j], limits[i]])
for ax in axs.flat:
    ax.set_xticks([])
    ax.set_yticks([])
plt.show()
