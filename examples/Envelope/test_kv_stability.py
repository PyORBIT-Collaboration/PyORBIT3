"""KV envelope stability versus bare and depressed phase advance in FODO.

Linearize the fixed-emittance Twiss map about each matched envelope. Its
four eigenvalues describe envelope modes, unlike the particle transfer map.
Scan kq and intensity, measure both phase advances, and interpolate the heatmap.
Blank regions lie outside the sampled domain; zero depressed phase is a limit.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from utils import make_envelope
from utils import make_fodo
from utils import make_output_dir
from utils import matched_envelope
from utils import twiss_parameters


def phase_advances(matrix: np.ndarray) -> np.ndarray:
    return np.degrees(
        [
            np.arccos(np.clip(np.trace(matrix[i : i + 2, i : i + 2]) / 2, -1, 1))
            for i in (0, 2)
        ]
    )


def main(args: argparse.Namespace) -> None:
    out = make_output_dir(__file__)
    plt.style.use("style.mplstyle")
    rows = []
    failures = 0
    for kq in np.linspace(args.kq_min, args.kq_max, args.points):
        lattice, _ = make_fodo(kq, args.nslice)
        try:
            bare, _ = matched_envelope(lattice, intensity=0.0, optimizer_kwargs={"verbose": 0})
        except (RuntimeError, ValueError):
            failures += args.points
            continue
        matrix = lattice.trackEnvelope(
            bare.copy(), fit=False, static=True, return_matrix=True
        )
        bare_phase = phase_advances(matrix)
        initial = twiss_parameters(bare)
        for intensity in np.linspace(0.0, args.intensity_max, args.points):
            try:
                beam, error = matched_envelope(
                    lattice, intensity, initial, optimizer_kwargs={"verbose": 0}
                )
            except (RuntimeError, ValueError):
                failures += 1
                continue
            initial = twiss_parameters(beam)

            def turn(parameters: np.ndarray) -> np.ndarray:
                trial = make_envelope(parameters, intensity)
                lattice.trackEnvelope(trial, fit=False, static=True, sc="2d")
                return twiss_parameters(trial)

            step = args.delta
            directions = step * np.eye(4)
            jacobian = np.column_stack(
                [
                    (turn(initial + d) - turn(initial - d)) / (2 * step)
                    for d in directions
                ]
            )
            growth = np.log(max(abs(np.linalg.eigvals(jacobian))))
            matrix = lattice.trackEnvelope(
                beam.copy(), fit=False, static=True, sc="2d", return_matrix=True
            )
            rows.append(
                [kq, intensity, *bare_phase, *phase_advances(matrix), growth, error]
            )
    data = np.asarray(rows)
    if not len(data):
        raise RuntimeError("No matched solutions found in scan")

    np.savetxt(
        out / "stability.csv",
        data,
        delimiter=",",
        header="kq,intensity,bare_x_deg,bare_y_deg,depressed_x_deg,depressed_y_deg,log_spectral_radius,match_residual",
    )
    print(f"Matched points: {len(data)}; skipped/failed points: {failures}")
    print(
        f"Unstable points (growth > {args.growth_tol:g} per cell): {sum(data[:, 6] > args.growth_tol)}"
    )
    phase_grid = np.linspace(0.0, 120.0, args.grid_points)
    bare_grid, depressed_grid = np.meshgrid(phase_grid, phase_grid)
    # The physics samples are uniform in kq and intensity, not phase advance.
    # Interpolate only for display; retain the actual samples in stability.csv.
    if np.linalg.matrix_rank(data[:, [2, 4]] - data[0, [2, 4]]) < 2:
        raise RuntimeError("Heatmap requires samples spanning both phase axes")
    growth_grid = griddata(
        data[:, [2, 4]], np.maximum(data[:, 6], 0.0),
        (bare_grid, depressed_grid), method="linear",
    )
    growth_grid[depressed_grid > bare_grid] = np.nan
    fig, ax = plt.subplots()
    image = ax.pcolormesh(
        bare_grid, depressed_grid, growth_grid, shading="auto", cmap="viridis"
    )
    ax.plot(phase_grid, phase_grid, color="white", linewidth=0.8, linestyle="--")
    ax.set(
        xlabel="Bare phase advance x [deg/cell]",
        ylabel="Depressed phase advance x [deg/cell]",
        xlim=(0.0, 120.0),
        ylim=(0.0, 120.0),
    )
    fig.colorbar(image, ax=ax, label="Envelope growth exponent [1/cell]")
    fig.savefig(out / "stability.png")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kq-min", type=float, default=0.1)
    parser.add_argument("--kq-max", type=float, default=1.15)
    parser.add_argument("--intensity-max", type=float, default=1e10)
    parser.add_argument("--points", type=int, default=31, help="Resolution of the kq-intensity scan")
    parser.add_argument("--grid-points", type=int, default=121, help="Interpolated heatmap resolution from 0 to 120 degrees")
    parser.add_argument("--nslice", type=int, default=20)
    parser.add_argument("--delta", type=float, default=1e-5)
    parser.add_argument("--growth-tol", type=float, default=1e-4)
    main(parser.parse_args())
