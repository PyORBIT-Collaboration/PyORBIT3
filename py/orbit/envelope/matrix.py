import math

import numpy as np

from ..core.bunch import Bunch
from ..core.bunch import SyncParticle
from ..lattice import AccNode
from ..teapot import DriftTEAPOT
from ..teapot import QuadTEAPOT
from ..teapot import BendTEAPOT
from ..teapot import TiltTEAPOT
from ..teapot import KickTEAPOT
from ..teapot import ApertureTEAPOT
from ..teapot import BunchWrapTEAPOT
from ..teapot import FringeFieldTEAPOT
from ..teapot import MonitorTEAPOT
from ..teapot import TurnCounterTEAPOT
from ..teapot import MultipoleTEAPOT
from ..teapot import NodeTEAPOT
from ..teapot import SolenoidTEAPOT
from ..utils import speed_of_light


def get_dp_p_coeff(sync_part: SyncParticle) -> float:
    beta = sync_part.beta()
    gamma = sync_part.gamma()
    rest_energy = sync_part.mass()  # GeV

    # dE/E = (beta^2) * dp/p
    # dE = (beta^2 * E) * dp/p
    # dE = (beta^2 * gamma * m * c^2) * dp/p
    return 1.0 / (beta**2 * gamma * rest_energy)


def convert_matrix_dp_p_to_dE(matrix: np.ndarray, sync_part: SyncParticle) -> np.ndarray:
    # v = [x, x', y, y', z, dp/p]
    # w = [x, x', y, y', z, dE]
    # v = A w
    # v -> M v
    # w -> A M A^-1

    # scale = np.identity(7)
    # scale[5, 5] = dp_p_coeff
    #
    # scale_inv = np.identity(7)
    # scale_inv[5, 5] = 1.0 / dp_p_coeff
    #
    # return np.linalg.multi_dot([scale, matrix, scale_inv])

    dp_p_coeff = get_dp_p_coeff(sync_part)
    matrix[:5, 5] *= dp_p_coeff
    matrix[5, :5] /= dp_p_coeff
    matrix[5, 6] /= dp_p_coeff
    return matrix


class MatrixFactory:
    """Factory for 7 x 7 transfer matrices.

    Units: x [m], x' [rad], y [m], y' [rad], z [m], dE [GeV]
    """
    def __init__(self, handle_unknown: str | None = None) -> None:
        self.ignore_node_types = [
            ApertureTEAPOT,
            BunchWrapTEAPOT,
            FringeFieldTEAPOT,
            MonitorTEAPOT,
            TurnCounterTEAPOT,
            NodeTEAPOT,
        ]
        self.handle_unknown = handle_unknown

    def drift_matrix(self, length: float, sync_part: SyncParticle) -> np.ndarray:
        M = np.identity(7)
        M[0, 1] = length
        M[2, 3] = length
        M[4, 5] = length / sync_part.gamma() ** 2
        M = convert_matrix_dp_p_to_dE(M, sync_part)
        return M

    def quad_matrix(self, length: float, kq: float, sync_part: SyncParticle) -> np.ndarray:
        if abs(kq) == 0:
            return self.drift_matrix(length=length, sync_part=sync_part)

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

        M[4, 5] = length / sync_part.gamma() ** 2
        M = convert_matrix_dp_p_to_dE(M, sync_part)
        return M

    def bend_matrix(self, length: float, theta: float, sync_part: SyncParticle) -> np.ndarray:
        if length <= 0:
            return np.identity(7)

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
        M = convert_matrix_dp_p_to_dE(M, sync_part)
        return M

    def tilt_matrix(self, angle: float) -> np.ndarray:
        M = np.identity(7)
        M[0, 0] = M[1, 1] = +math.cos(angle)
        M[0, 2] = M[1, 3] = -math.sin(angle)
        M[2, 0] = M[3, 1] = +math.sin(angle)
        M[2, 2] = M[3, 3] = +math.cos(angle)
        return M

    def translation_matrix(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> np.ndarray:
        M = np.identity(7)
        M[0, -1] = x
        M[2, -1] = y
        M[4, -1] = z
        return M

    def kick_matrix(self, kx: float = 0.0, ky: float = 0.0, kE: float = 0.0) -> np.ndarray:
        M = np.identity(7)
        M[1, -1] = kx
        M[3, -1] = ky
        M[5, -1] = kE
        return M

    def solenoid_matrix(self, length: float, B: float, sync_part: SyncParticle) -> np.ndarray:
        if B == 0:
            return self.drift_matrix(length=length, sync_part=sync_part)

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
        M[3, 2] = math.sin(phase) * (-B)
        M[3, 3] = math.cos(phase)
        M[4, 5] = length / sync_part.gamma() ** 2

        M = np.linalg.multi_dot([np.linalg.inv(V), M, V])
        M = convert_matrix_dp_p_to_dE(M, sync_part)
        return M

    def cf_matrix(self, length: float, kq: float, sync_part: SyncParticle) -> np.ndarray:
        raise NotImplementedError()

    def __call__(self, node: AccNode, sync_part: SyncParticle, part_index: int = 0) -> np.ndarray:
        if type(node) is DriftTEAPOT:
            length = node.getLength(part_index)
            return self.drift_matrix(length=length, sync_part=sync_part)

        elif type(node) is QuadTEAPOT:
            length = node.getLength(part_index)

            scale = 1.0
            if node.waveform:
                scale = node.waveform.getStrength()

            kq = scale * node.getParam("kq")
            return self.quad_matrix(length=length, kq=kq, sync_part=sync_part)

        elif type(node) is BendTEAPOT:
            length = node.getLength(part_index)
            theta = node.getParam("theta") / node.getnParts()
            return self.bend_matrix(length=length, theta=theta, sync_part=sync_part)

        elif type(node) is KickTEAPOT:
            length = node.getLength(part_index)
            nparts = node.getnParts()

            scale = 1.0
            if node.waveform is not None:
                scale = node.waveform.getStrength()

            kx = scale * node.getParam("kx") / nparts
            ky = scale * node.getParam("ky") / nparts
            kE = node.getParam("dE") / nparts

            return np.matmul(
                self.kick_matrix(kx=kx, ky=ky, kE=kE),
                self.drift_matrix(length=length, sync_part=sync_part)
            )

        elif type(node) is TiltTEAPOT:
            angle = node.getTiltAngle()
            return self.tilt_matrix(angle)

        elif type(node) is SolenoidTEAPOT:
            B = node.getParam("B")
            if node.waveform is not None:
                B *= node.waveform.getStrength()
            length = node.getLength(part_index)
            return self.solenoid_matrix(length=length, B=B, sync_part=sync_part)

        elif type(node) in self.ignore_node_types:
            return np.identity(7)

        else:
            if type(node) is MultipoleTEAPOT:
                if np.all(np.abs(node.getParam("kls")) == 0):
                    return self.drift_matrix(length=node.getLength(), sync_part=sync_part)

            elif self.handle_unknown == "drift":
                return self.drift_matrix(length=node.getLength(), sync_part=sync_part)

            elif self.handle_unknown == "fit":
                raise NotImplementedError()

            raise NotImplementedError("Unsupported node: {}.".format(node))
