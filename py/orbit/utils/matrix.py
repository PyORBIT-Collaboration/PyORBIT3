"""Transfer matrix definitions."""
import math

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.core.orbit_utils import Matrix
from orbit.core.teapot_base import MatrixGenerator
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import ApertureTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import BendTEAPOT
from orbit.teapot import KickTEAPOT
from orbit.teapot import MonitorTEAPOT
from orbit.teapot import MultipoleTEAPOT
from orbit.teapot import NodeTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.teapot import SolenoidTEAPOT
from orbit.teapot import FringeFieldTEAPOT
from orbit.teapot import BunchWrapTEAPOT
from orbit.teapot import TiltTEAPOT
from orbit.teapot import ContinuousLinearFocusingTEAPOT
from orbit.teapot import TurnCounterTEAPOT
from orbit.py_linac.lattice import MarkerLinacNode as MarkerLINAC
from orbit.py_linac.lattice import Drift as DriftLINAC
from orbit.py_linac.lattice import Quad as QuadLINAC
from orbit.py_linac.lattice import Bend as BendLINAC
from orbit.py_linac.lattice import DCorrectorH as DCorrectorHLINAC
from orbit.py_linac.lattice import DCorrectorV as DCorrectorVLINAC
from orbit.py_linac.lattice import Solenoid as SolenoidLINAC
from orbit.py_linac.lattice import TiltElement as TiltLINAC
from orbit.py_linac.lattice import FringeField as FringeFieldLINAC
from orbit.py_linac.lattice import BaseRF_Gap as BaseRF_Gap
from orbit.utils.consts import speed_of_light


IGNORE_NODE_TYPES = [
    NodeTEAPOT,
    MonitorTEAPOT,
    FringeFieldTEAPOT,
    ApertureTEAPOT,
    BunchWrapTEAPOT,
    TurnCounterTEAPOT,
    MarkerLINAC,
    FringeFieldLINAC,
]


def get_dp_p_coeff(sync_part: SyncParticle) -> float:
    # dE/E = (beta^2) * dp/p
    # dE = (beta^2 * E) * dp/p
    # dE = (beta^2 * gamma * m * c^2) * dp/p
    beta = sync_part.beta()
    gamma = sync_part.gamma()
    rest_energy = sync_part.mass()  # GeV
    return 1.0 / (beta**2 * gamma * rest_energy)


def get_zp_coeff(sync_part: SyncParticle) -> float:
    # dE/E = (beta^2) * dp/p = (beta^2) * (gamma^2) z'
    # dE = (beta^2 * gamma^2 * E) * z'
    # dE = (beta^2 * gamma^3 * m * c^2) * z'
    beta = sync_part.beta()
    gamma = sync_part.gamma()
    rest_energy = sync_part.mass()
    return 1.0 / (beta**2 * gamma**3 * rest_energy)


def convert_matrix_dp_p_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    # v = [x, x', y, y', z, dp/p]
    # w = [x, x', y, y', z, dE]
    # v = A w
    # v -> M v
    # w -> A M A^-1
    dp_p_coeff = get_dp_p_coeff(sync_part)
    matrix[:5, 5] *= dp_p_coeff
    matrix[5, :5] /= dp_p_coeff
    matrix[5, 6] /= dp_p_coeff  # driving term
    return matrix


def convert_matrix_zp_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    zp_coeff = get_zp_coeff(sync_part)
    matrix[:5, 5] *= zp_coeff
    matrix[5, :5] /= zp_coeff
    matrix[5, 6] /= zp_coeff  # driving term
    return matrix


def get_matrix_tilt(angle: float) -> np.ndarray:
    cos_phi = math.cos(angle)
    sin_phi = math.sin(angle)

    M = np.identity(7)
    M[0, 0] = M[1, 1] = +cos_phi
    M[0, 2] = M[1, 3] = -sin_phi
    M[2, 0] = M[3, 1] = +sin_phi
    M[2, 2] = M[3, 3] = +cos_phi
    return M


def get_matrix_kick(kx: float = 0.0, ky: float = 0.0, kE: float = 0.0) -> np.ndarray:
    M = np.identity(7)
    M[1, -1] = kx
    M[3, -1] = ky
    M[5, -1] = kE
    return M


def get_matrix_drift(sync_part: SyncParticle, length: float) -> np.ndarray:
    M = np.identity(7)
    M[0, 1] = length
    M[2, 3] = length
    M[4, 5] = length / (sync_part.gamma() ** 2)
    M[4, 5] *= get_dp_p_coeff(sync_part)  # convert_matrix_dp_p_to_dE(M, sync_part)

    sync_part.time(sync_part.time() + length / (sync_part.beta() * speed_of_light))
    return M


def get_matrix_quad(sync_part: SyncParticle, length: float, kq: float) -> np.ndarray:
    if abs(kq) == 0:
        return get_matrix_drift(sync_part, length)

    sqrt_abs_kq = math.sqrt(abs(kq))

    M = np.identity(7)
    if kq > 0:
        cx = np.cos(sqrt_abs_kq * length)
        sx = np.sin(sqrt_abs_kq * length)
        cy = np.cosh(sqrt_abs_kq * length)
        sy = np.sinh(sqrt_abs_kq * length)
        M[0, 0] = cx
        M[0, 1] = +sx / sqrt_abs_kq
        M[1, 0] = -sx * sqrt_abs_kq
        M[1, 1] = cx
        M[2, 2] = cy
        M[2, 3] = sy / sqrt_abs_kq
        M[3, 2] = sy * sqrt_abs_kq
        M[3, 3] = cy
    elif kq < 0:
        cx = np.cosh(sqrt_abs_kq * length)
        sx = np.sinh(sqrt_abs_kq * length)
        cy = np.cos(sqrt_abs_kq * length)
        sy = np.sin(sqrt_abs_kq * length)
        M[0, 0] = cx
        M[0, 1] = sx / sqrt_abs_kq
        M[1, 0] = sx * sqrt_abs_kq
        M[1, 1] = cx
        M[2, 2] = cy
        M[2, 3] = +sy / sqrt_abs_kq
        M[3, 2] = -sy * sqrt_abs_kq
        M[3, 3] = cy

    M[4, 5] = length / (sync_part.gamma()**2)
    M[4, 5] *= get_dp_p_coeff(sync_part)  # convert_matrix_dp_p_to_dE(M, sync_part)

    sync_part.time(sync_part.time() + length / (sync_part.beta() * speed_of_light))
    return M


def get_matrix_bend(sync_part: SyncParticle, length: float, theta: float) -> np.ndarray:
    rho = length / theta
    cx = math.cos(theta)
    sx = math.sin(theta)

    M = np.identity(7)
    M[0, 0] = cx
    M[0, 1] = rho * sx
    M[0, 5] = rho * (1.0 - cx)
    M[1, 0] = -sx / rho
    M[1, 1] = cx
    M[1, 5] = sx
    M[2, 3] = length
    M[4, 0] = -sx
    M[4, 1] = -rho * (1.0 - cx)
    M[4, 5] = -(sync_part.beta() ** 2) * length + rho * sx
    M[:5, 5] *= get_dp_p_coeff(sync_part)  # convert_matrix_dp_p_to_dE(M, sync_part)

    sync_part.time(sync_part.time() + length / (sync_part.beta() * speed_of_light))
    return M


def get_matrix_solenoid(sync_part: SyncParticle, length: float, B: float) -> np.ndarray:
    if B == 0:
        return get_matrix_drift(sync_part, length)

    phase = B * length

    V = np.identity(7)
    V[:4, :4] = 0.0
    V[0, 1] = -1.0 / B
    V[0, 2] = 0.5
    V[1, 0] = 0.5 * B
    V[1, 3] = 1.0
    V[2, 1] = 1.0 / B
    V[2, 2] = 0.5
    V[3, 0] = -0.5 * B
    V[3, 3] = 1.0

    M = np.identity(7)
    M[0, 0] = +1.0
    M[1, 1] = -1.0
    M[2, 2] = math.cos(phase)
    M[2, 3] = math.sin(phase) / B
    M[3, 2] = math.sin(phase) * B * -1.0
    M[3, 3] = math.cos(phase)
    M[4, 5] = length / (sync_part.gamma()**2)

    M = np.linalg.inv(V) @ M @ V
    M[4, 5] *= get_dp_p_coeff(sync_part)  # convert_matrix_dp_p_to_dE(M, sync_part)

    sync_part.time(sync_part.time() + length / (sync_part.beta() * speed_of_light))
    return M


def get_matrix_cf(sync_part: SyncParticle, length: float, kq: float) -> np.ndarray:
    if kq == 0:
        return get_matrix_drift(sync_part, length)

    sqrt_abs_kq = math.sqrt(abs(kq))

    cx = math.cos(sqrt_abs_kq * length)
    sx = math.sin(sqrt_abs_kq * length)

    M = np.identity(7)
    M[0, 0] = M[2, 2] = cx
    M[0, 1] = M[2, 3] = +sx / sqrt_abs_kq
    M[1, 0] = M[3, 2] = -sx * sqrt_abs_kq
    M[1, 1] = M[3, 3] = cx
    M[4, 5] = length / (sync_part.gamma()**2)
    M[4, 5] *= get_dp_p_coeff(sync_part)

    sync_part.time(sync_part.time() + length / (sync_part.beta() * speed_of_light))
    return M


def get_matrix_rf_gap(sync_part: SyncParticle, frequency: float, E0TL: float, phase: float) -> np.ndarray:
    gamma = sync_part.gamma()
    beta = sync_part.beta()
    mass = sync_part.mass()
    charge = sync_part.charge()

    kin_energy_in = sync_part.kinEnergy()
    charge_E0TL_sin = charge * E0TL * math.sin(phase)
    kin_energy_delta = charge * E0TL * math.cos(phase)

    # Calculate parameters in the center of the gap.
    sync_part.momentum(sync_part.energyToMomentum(kin_energy_in + kin_energy_delta / 2.0))
    gamma_gap = sync_part.gamma()
    beta_gap = sync_part.beta()

    # Move to the end of the gap.
    kin_energy_out = kin_energy_in + kin_energy_delta
    sync_part.momentum(sync_part.energyToMomentum(kin_energy_out))

    # The base RF gap is simple - no phase correction.
    gamma_out = sync_part.gamma()
    beta_out = sync_part.beta()
    prime_coeff = (beta * gamma) / (beta_out * gamma_out)

    # Wave momentum
    k = 2.0 * math.pi * frequency / speed_of_light
    phase_time_coeff = k / beta

    # Transverse focusing coefficient
    kappa = -charge * E0TL * k / (2.0 * mass * beta_gap**2 * beta_out * gamma_gap**2 * gamma_out)
    d_rp = kappa * math.sin(phase)

    M = np.eye(7)
    M[5, 4] = charge_E0TL_sin * phase_time_coeff
    M[4, 4] = beta_out / beta
    M[1, 1] = prime_coeff
    M[3, 3] = prime_coeff
    M[1, 0] = d_rp
    M[3, 2] = d_rp
    return M


def orbit_matrix_to_numpy(matrix: Matrix) -> np.ndarray:
    matrix_out = np.zeros(matrix.size())
    for i in range(matrix_out.shape[0]):
        for j in range(matrix_out.shape[1]):
            matrix_out[i, j] = matrix.get(i, j)
    return matrix_out


def get_matrix(node: AccNode, sync_part: SyncParticle, part_index: int = -1, fit: bool = False) -> np.ndarray | None:
    """Calculate transfer matrix and update synchronous particle.

    This function maps various accelerator nodes to 7 x 7 transfer matrices.
    For non-accelerating, finite-length nodes, the synchronous particle time
    is updated as in a drift. Accelerating nodes such as RF gaps will update
    the synchronous particle energy.

    Args:
        node: The accelerator node.
        sync_part: Synchronous particle
        part_index: Index of the part within the node. An index of -1 returns
            the transfer matrix for the entire node.
        fit: Whether to fit matrix to tracking data.
    Returns:
        7 x 7 transfer matrix or None.
    """
    if fit:
        raise NotImplementedError

    node_type = type(node)
    if node_type in IGNORE_NODE_TYPES:
        return None

    length = node.getLength(part_index)
    nparts = node.getnParts()

    if node_type is DriftTEAPOT:
        if length <= 0:
            return None
        return get_matrix_drift(sync_part, length)

    elif node_type is SolenoidTEAPOT:
        if length <= 0:
            return None

        B = node.getParam("B")
        if node.waveform:
            B *= node.waveform.getStrength()
        B *= np.sign(sync_part.charge())

        return get_matrix_solenoid(sync_part, length=length, B=B)

    elif node_type is MultipoleTEAPOT:
        if length <= 0:
            return None

        if np.all(np.abs(node.getParam("kls")) == 0):
            return get_matrix_drift(sync_part, length)

    elif node_type is QuadTEAPOT:
        if length <= 0:
            return None

        kq = node.getParam("kq")
        if node.waveform:
            kq *= node.waveform.getStrength()
        kq *= np.sign(sync_part.charge())

        return get_matrix_quad(sync_part, length=length, kq=kq)

    elif node_type is BendTEAPOT:
        if length <= 0:
            return None

        theta = node.getParam("theta") / (nparts - 1)
        if part_index == 0 or part_index == nparts - 1:
            theta *= 0.5
        theta *= np.sign(sync_part.charge())

        return get_matrix_bend(sync_part, length=length, theta=theta)

    elif node_type is KickTEAPOT:
        scale = 1.0
        if node.waveform is not None:
            scale = node.waveform.getStrength()

        scale /= (nparts - 1)
        kx = scale * node.getParam("kx")
        ky = scale * node.getParam("ky")
        kE = node.getParam("dE")

        if abs(kx) > 0 or abs(ky) > 0 or abs(kE) > 0:
            return np.matmul(
                get_matrix_kick(kx=kx, ky=ky, kE=kE),
                get_matrix_drift(sync_part, length),
            )
        else:
            return get_matrix_drift(sync_part, length)

    elif node_type is TiltTEAPOT:
        angle = node.getTiltAngle()
        if angle == 0:
            return None
        return get_matrix_tilt(angle)

    elif node_type is ContinuousLinearFocusingTEAPOT:
        if length <= 0:
            return None

        kq = node.getParam("kq")
        kq *= np.sign(sync_part.charge())
        if node.waveform:
            kq *= node.waveform.getStrength()

        return get_matrix_cf(sync_part, length=length, kq=kq)

    elif node_type is DriftLINAC:
        if length <= 0:
            return None
        return get_matrix_drift(sync_part, length)

    elif node_type is QuadLINAC:
        if length <= 0:
            return None

        brho = 3.335640952 * sync_part.momentum() / sync_part.charge()
        kq = node.getParam("dB/dr") / brho
        return get_matrix_quad(sync_part, length=length, kq=kq)

    elif node_type is BendLINAC:
        if length <= 0:
            return None

        theta = node.getParam("theta") / (nparts - 1)
        if part_index == 0 or part_index == nparts - 1:
            theta *= 0.5
        theta *= np.sign(sync_part.charge())

        return get_matrix_bend(sync_part, length=length, theta=theta)

    elif node_type is DCorrectorHLINAC:
        length = node.getParam("effLength") / nparts
        field = node.getParam("B")
        delta_xp = -field * sync_part.charge() * length * 0.299792 / sync_part.momentum()
        if delta_xp == 0:
            return None
        return get_matrix_kick(kx=delta_xp, ky=0.0, kE=0.0)

    elif node_type is DCorrectorVLINAC:
        length = node.getParam("effLength") / nparts
        field = node.getParam("B")
        delta_yp = -field * sync_part.charge() * length * 0.299792 / sync_part.momentum()
        if delta_yp == 0:
            return None
        return get_matrix_kick(kx=0.0, ky=delta_yp, kE=0.0)

    elif node_type is SolenoidLINAC:
        if length <= 0:
            return None
        B = node.getParam("B") * np.sign(sync_part.charge())
        return get_matrix_solenoid(sync_part, length=length, B=B)

    elif node_type is TiltLINAC:
        angle = node.getTiltAngle()
        if angle == 0:
            return None
        return get_matrix_tilt(angle=angle)

    elif node_type is BaseRF_Gap:
        E0TL = node.getParam("E0TL")
        mode_phase = node.getParam("mode") * math.pi

        cavity = node.getRF_Cavity()
        frequency = cavity.getFrequency()
        phase = cavity.getPhase() + mode_phase
        amplitude = cavity.getAmp()

        arrival_time = sync_part.time()
        arrival_time_design = cavity.getDesignArrivalTime()

        if node.isFirstRFGap():
            if cavity.isDesignSetUp():
                phase = math.fmod(frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase, 2.0 * math.pi)
            else:
                raise ValueError("Run `trackDesign` first to initialize cavity phases.")
        else:
            phase = math.fmod(frequency * (arrival_time - arrival_time_design) * 2.0 * math.pi + phase,2.0 * math.pi)

        node.setGapPhase(phase)

        if amplitude == 0.0:
            return None

        return get_matrix_rf_gap(
            sync_part=sync_part,
            frequency=frequency,
            E0TL=(E0TL * amplitude),
            phase=phase,
        )

    raise NotImplementedError(str(node))