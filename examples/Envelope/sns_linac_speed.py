"""Profile envelope tracking in SNS linac."""

import argparse
import cProfile
import pstats

from orbit.bunch_generators import GaussDist3D
from orbit.bunch_generators import KVDist3D
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D
from orbit.core.bunch import Bunch
from orbit.core.linac import MatrixRfGap
from orbit.envelope import Envelope
from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.utils.consts import charge_electron
from orbit.utils.consts import mass_electron
from orbit.utils.consts import mass_proton


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist", choices=("kv", "waterbag", "gauss"), default="kv")
    parser.add_argument("--nparts", type=int, default=20_000)
    parser.add_argument("--current", type=float, default=0.038)
    parser.add_argument("--turns", type=int, default=5)
    parser.add_argument("--seq-stop", default="CCL2")
    parser.add_argument("--sc", type=int, default=0)
    parser.add_argument("--sc-model", choices=("ellipsoid", "fft"), default="ellipsoid")
    parser.add_argument("--sc-grid", type=int, default=64)
    parser.add_argument("--sc-ellipsoid-n", type=int, default=1)
    parser.add_argument("--sc-path-length-min", type=float, default=0.01)
    parser.add_argument("--fit", type=int, default=0)
    parser.add_argument("--fallback", type=int, default=0)
    parser.add_argument("--static", type=int, default=0)
    parser.add_argument(
        "--bunch", action="store_true", help="Also benchmark particle tracking"
    )
    return parser.parse_args()


def make_lattice(args: argparse.Namespace) -> LinacAccLattice:
    sequences = [
        "MEBT",
        "DTL1",
        "DTL2",
        "DTL3",
        "DTL4",
        "DTL5",
        "DTL6",
        "CCL1",
        "CCL2",
        "CCL3",
        "CCL4",
        "SCLMed",
        "SCLHigh",
        "HEBT1",
        "HEBT2",
    ]
    if args.seq_stop:
        sequences = sequences[: sequences.index(args.seq_stop) + 1]

    factory = SNS_LinacLatticeFactory()
    factory.setMaxDriftLength(args.sc_path_length_min)
    lattice = factory.getLinacAccLattice(sequences, "inputs/sns_linac.xml")
    for node in lattice.getNodes():
        try:
            node.setUsageFringeFieldIN(False)
            node.setUsageFringeFieldOUT(False)
        except AttributeError:
            pass
    for gap in lattice.getRF_Gaps():
        gap.setCppGapModel(MatrixRfGap())
    return lattice


def make_bunch(args: argparse.Namespace) -> Bunch:
    kinetic_energy = 0.0025
    mass = mass_proton + 2 * mass_electron
    charge = -1.0
    intensity = args.current / 402.5e6 / (abs(charge) * charge_electron)

    bunch = Bunch()
    bunch.mass(mass)
    bunch.charge(charge)
    bunch.macroSize(intensity / args.nparts)
    bunch.getSyncParticle().kinEnergy(kinetic_energy)
    bunch.getSyncParticle().time(0.0)

    twiss = [
        TwissContainer(-1.962, 0.183, 2.874e-06),
        TwissContainer(+1.768, 0.162, 2.874e-06),
        TwissContainer(-0.0196, 116.414, 1.651e-08),
    ]
    if args.dist == "kv":
        distribution = KVDist3D(*twiss)
    elif args.dist == "waterbag":
        distribution = WaterBagDist3D(*twiss)
    elif args.dist == "gauss":
        distribution = GaussDist3D(*twiss)
    else:
        raise ValueError
    for _ in range(args.nparts):
        bunch.addParticle(*distribution.getCoordinates())
    return bunch


def main(args: argparse.Namespace) -> None:
    lattice = make_lattice(args)
    bunch = make_bunch(args)
    lattice.trackDesignBunch(bunch)
    envelope = Envelope(bunch=bunch)

    profiler = cProfile.Profile()
    profiler.enable()
    lattice.trackEnvelope(
        envelope,
        sc=("3d" if args.sc else None),
        fit=args.fit,
        static=args.static,
        fallback=args.fallback,
    )
    profiler.disable()

    profiler_stats = pstats.Stats(profiler)
    profiler_stats.sort_stats(pstats.SortKey.TIME)
    profiler_stats.print_stats(10)


if __name__ == "__main__":
    main(parse_args())
