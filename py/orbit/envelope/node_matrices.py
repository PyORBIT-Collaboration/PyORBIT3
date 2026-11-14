import math
from functools import singledispatch

import numpy as np

from orbit.teapot.teapot import TurnCounterTEAPOT
from orbit.teapot.teapot import NodeTEAPOT
from orbit.teapot.teapot import DriftTEAPOT
from orbit.teapot.teapot import ApertureTEAPOT
from orbit.teapot.teapot import MonitorTEAPOT
from orbit.teapot.teapot import BunchWrapTEAPOT
from orbit.teapot.teapot import SolenoidTEAPOT
from orbit.teapot.teapot import MultipoleTEAPOT
from orbit.teapot.teapot import QuadTEAPOT
from orbit.teapot.teapot import BendTEAPOT
from orbit.teapot.teapot import KickTEAPOT
from orbit.teapot.teapot import TiltTEAPOT
from orbit.teapot.teapot import FringeFieldTEAPOT
from orbit.teapot.teapot import ContinuousLinearFocusingTEAPOT
from orbit.py_linac.lattice.LinacAccNodes import MarkerLinacNode
from orbit.py_linac.lattice.LinacAccNodes import Drift
from orbit.py_linac.lattice.LinacAccNodes import Quad
from orbit.py_linac.lattice.LinacAccNodes import Bend
from orbit.py_linac.lattice.LinacAccNodes import DCorrectorH
from orbit.py_linac.lattice.LinacAccNodes import DCorrectorV
from orbit.py_linac.lattice.LinacAccNodes import Solenoid
from orbit.py_linac.lattice.LinacAccNodes import TiltElement
from orbit.py_linac.lattice.LinacAccNodes import FringeField
from orbit.py_linac.lattice.LinacRfGapNodes import BaseRF_Gap
from orbit.space_charge.sc3d import SC3D_AccNode
from orbit.space_charge.sc3d import SC_UniformEllipses_AccNode
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2dslicebyslice import SC2DSliceBySlice_AccNode

from .envelope import Envelope
from .matrix import get_matrix_bend
from .matrix import get_matrix_bend_edge
from .matrix import get_matrix_cf
from .matrix import get_matrix_drift
from .matrix import get_matrix_kick
from .matrix import get_matrix_quad
from .matrix import get_matrix_rf_gap
from .matrix import get_matrix_solenoid
from .matrix import get_matrix_tilt


@singledispatch
def get_node_matrix(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    raise NotImplementedError(str(node))


@get_node_matrix.register(ApertureTEAPOT)
@get_node_matrix.register(BunchWrapTEAPOT)
@get_node_matrix.register(MonitorTEAPOT)
@get_node_matrix.register(TurnCounterTEAPOT)
@get_node_matrix.register(MarkerLinacNode)
@get_node_matrix.register(SC3D_AccNode)
@get_node_matrix.register(SC_UniformEllipses_AccNode)
@get_node_matrix.register(SC2p5D_AccNode)
@get_node_matrix.register(SC2DSliceBySlice_AccNode)
def _matrix_identity(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    return None


@get_node_matrix.register(NodeTEAPOT)
def _matrix_node_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    if type(node) is NodeTEAPOT:
        return None
    raise NotImplementedError(str(node))


@get_node_matrix.register(DriftTEAPOT)
@get_node_matrix.register(Drift)
def _matrix_drift(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    return get_matrix_drift(envelope, length)


@get_node_matrix.register(SolenoidTEAPOT)
def _matrix_solenoid_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    B = node.getParam("B")
    if node.waveform:
        B *= node.waveform.getStrength()
    B *= np.sign(envelope.charge)
    return get_matrix_solenoid(envelope, length=length, B=B)


@get_node_matrix.register(MultipoleTEAPOT)
def _matrix_multipole_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    if np.all(np.abs(node.getParam("kls")) == 0):
        return get_matrix_drift(envelope, length)
    raise NotImplementedError(str(node))


@get_node_matrix.register(QuadTEAPOT)
def _matrix_quad_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    kq = node.getParam("kq")
    if node.waveform:
        kq *= node.waveform.getStrength()
    return get_matrix_quad(envelope, length=length, kq=kq)


@get_node_matrix.register(BendTEAPOT)
def _matrix_bend_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    if np.any(np.asarray(node.getParam("kls")) != 0):
        raise NotImplementedError

    length = node.getLength(part_index)
    if length <= 0:
        return None

    nparts = node.getnParts()
    if part_index < 0:
        theta = node.getParam("theta")
    else:
        theta = node.getParam("theta") / (nparts - 1)
        if part_index == 0 or part_index == nparts - 1:
            theta *= 0.5
    return get_matrix_bend(envelope, length=length, theta=theta)


@get_node_matrix.register(KickTEAPOT)
def _matrix_kick_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    nparts = node.getnParts()

    if part_index == nparts - 1:
        return get_matrix_drift(envelope, length) if length > 0 else None

    scale = 1.0
    if node.waveform is not None:
        scale = node.waveform.getStrength()

    kx = scale * node.getParam("kx") / (nparts - 1)
    ky = scale * node.getParam("ky") / (nparts - 1)
    kE = node.getParam("dE") / (nparts - 1)

    if abs(kx) == abs(ky) == abs(kE) == 0:
        return get_matrix_drift(envelope, length) if length > 0 else None

    matrix_kick = get_matrix_kick(kx=kx, ky=ky, kE=kE)
    matrix_drift = get_matrix_drift(envelope, length)
    return matrix_kick @ matrix_drift


@get_node_matrix.register(TiltTEAPOT)
@get_node_matrix.register(TiltElement)
def _matrix_tilt(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    angle = node.getTiltAngle()
    if angle == 0:
        return None
    return get_matrix_tilt(angle)


@get_node_matrix.register(FringeFieldTEAPOT)
@get_node_matrix.register(FringeField)
def _matrix_fringe_field(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    if node.getType() in ("bend fringe in", "bend fringe out"):
        usage = node.getUsage()
        # Linac fringes skip their entire callback when disabled; TEAPOT
        # fringes still apply the geometric wedge transformations.
        if isinstance(node, FringeField) and not usage:
            return None
        if np.any(np.asarray(node.getParam("kls")) != 0):
            raise NotImplementedError
        angle = node.getParam("ea1" if node.getType() == "bend fringe in" else "ea2")
        if angle == 0.0:
            return None
        return get_matrix_bend_edge(angle, node.getParam("rho"), fringe=usage)
    params_dict = node.getParamsDict()
    if params_dict.get("ea1", 0.0) != 0.0 or params_dict.get("ea2", 0.0) != 0.0:
        raise NotImplementedError("Bend edge fringes require a fitted matrix")
    return None


@get_node_matrix.register(ContinuousLinearFocusingTEAPOT)
def _matrix_continuous_focusing_teapot(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    kq = node.getParam("kq")
    if node.waveform:
        kq *= node.waveform.getStrength()
    return get_matrix_cf(envelope, length=length, kq=kq)


@get_node_matrix.register(Quad)
def _matrix_quad_linac(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    kq = node.getParam("dB/dr") / envelope.rigidity
    return get_matrix_quad(envelope, length=length, kq=kq)


@get_node_matrix.register(Bend)
def _matrix_bend_linac(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    if np.any(np.asarray(node.getParam("kls")) != 0):
        raise NotImplementedError("Combined-function bends require a fitted matrix")
    length = node.getLength(part_index)
    if length <= 0:
        return None
    nparts = node.getnParts()
    theta = node.getParam("theta") / (nparts - 1)
    if part_index == 0 or part_index == nparts - 1:
        theta *= 0.5
    return get_matrix_bend(envelope, length=length, theta=theta)


@get_node_matrix.register(DCorrectorH)
def _matrix_corrector_h(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getParam("effLength") / node.getnParts()
    field = node.getParam("B")
    delta_xp = -field * envelope.charge * length * 0.299792 / envelope.momentum
    if delta_xp == 0:
        return None
    return get_matrix_kick(kx=delta_xp, ky=0.0, kE=0.0)


@get_node_matrix.register(DCorrectorV)
def _matrix_corrector_v(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getParam("effLength") / node.getnParts()
    field = node.getParam("B")
    delta_yp = -field * envelope.charge * length * 0.299792 / envelope.momentum
    if delta_yp == 0:
        return None
    return get_matrix_kick(kx=0.0, ky=delta_yp, kE=0.0)


@get_node_matrix.register(Solenoid)
def _matrix_solenoid_linac(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    length = node.getLength(part_index)
    if length <= 0:
        return None
    B = node.getParam("B") * np.sign(envelope.charge)
    return get_matrix_solenoid(envelope, length=length, B=B)


@get_node_matrix.register(BaseRF_Gap)
def _matrix_rf_gap(node, envelope: Envelope, part_index: int = -1) -> np.ndarray | None:
    E0TL = node.getParam("E0TL")
    mode_phase = node.getParam("mode") * math.pi

    cavity = node.getRF_Cavity()
    frequency = cavity.getFrequency()
    phase = cavity.getPhase() + mode_phase
    amplitude = cavity.getAmp()

    arrival_time = envelope.time
    arrival_time_design = cavity.getDesignArrivalTime()

    if node.isFirstRFGap():
        if cavity.isDesignSetUp():
            phase = math.fmod(
                frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase,
                2.0 * math.pi,
            )
        else:
            raise ValueError("Run `trackDesign` first to initialize cavity phases.")
    else:
        phase = math.fmod(
            frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase,
            2.0 * math.pi,
        )

    node.setGapPhase(phase)

    if amplitude == 0.0:
        return None

    return get_matrix_rf_gap(
        envelope=envelope,
        frequency=frequency,
        E0TL=(E0TL * amplitude),
        phase=phase,
    )
