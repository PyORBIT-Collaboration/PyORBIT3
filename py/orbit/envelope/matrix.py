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

from ..utils import speed_of_light


def get_dp_p_coeff(sync_part: SyncParticle) -> float:
    return 1.0 / (sync_part.momentum() * sync_part.beta())


class MatrixFactory:
    """Factory for 7 x 7 transfer matrices.

    Units: x [m], x' [rad], y [m], y' [rad], z [m], dE [GeV]
    """
    def __init__(self, ignore_unknown: bool = False) -> None:
        self.ignore_node_types = [
            ApertureTEAPOT,
            BunchWrapTEAPOT,
            FringeFieldTEAPOT,
            MonitorTEAPOT,
            TurnCounterTEAPOT,
        ]
        self.ignore_unknown = ignore_unknown

    def drift(self, length: float, sync_part: SyncParticle) -> np.ndarray:
        matrix = np.identity(7)
        matrix[0, 1] = length
        matrix[2, 3] = length
        matrix[4, 5] = length / sync_part.gamma() ** 2
        matrix[4, 5] *= get_dp_p_coeff(sync_part)
        return matrix

    def quad(self, length: float, kq: float, sync_part: SyncParticle) -> np.ndarray:
        sqrt_abs_kq = math.sqrt(abs(kq))

        matrix = np.identity(7)
        if kq > 0:
            cx = np.cos(sqrt_abs_kq * length)
            sx = np.sin(sqrt_abs_kq * length)
            cy = np.cosh(sqrt_abs_kq * length)
            sy = np.sinh(sqrt_abs_kq * length)
            matrix[0, 0] = cx
            matrix[0, 1] = +sx / sqrt_abs_kq
            matrix[1, 0] = -sx * sqrt_abs_kq
            matrix[1, 1] = cx
            matrix[2, 2] = cy
            matrix[2, 3] = sy / sqrt_abs_kq
            matrix[3, 2] = sy * sqrt_abs_kq
            matrix[3, 3] = cy
        elif kq < 0:
            cx = np.cosh(sqrt_abs_kq * length)
            sx = np.sinh(sqrt_abs_kq * length)
            cy = np.cos(sqrt_abs_kq * length)
            sy = np.sin(sqrt_abs_kq * length)
            matrix[0, 0] = cx
            matrix[0, 1] = sx / sqrt_abs_kq
            matrix[1, 0] = sx * sqrt_abs_kq
            matrix[1, 1] = cx
            matrix[2, 2] = cy
            matrix[2, 3] = +sy / sqrt_abs_kq
            matrix[3, 2] = -sy * sqrt_abs_kq
            matrix[3, 3] = cy

        matrix[4, 5] = length / sync_part.gamma() ** 2
        matrix[4, 5] *= get_dp_p_coeff(sync_part)
        return matrix

    def bend(self, length: float, theta: float, sync_part: SyncParticle) -> np.ndarray:
        if length <= 0:
            return np.identity(7)

        v = speed_of_light * sync_part.beta()
        sync_part.time(sync_part.time() + length / v)

        rho = length / theta
        cx = math.cos(theta)
        sx = math.sin(theta)

        matrix = np.identity(7)
        matrix[0, 0] = cx
        matrix[0, 1] = rho * sx
        matrix[0, 5] = rho * (1.0 - cx)
        matrix[1, 0] = -sx / rho
        matrix[1, 1] = cx
        matrix[1, 5] = sx

        matrix[2, 3] = length

        matrix[4, 0] = -sx
        matrix[4, 1] = -rho * (1.0 - cx)
        matrix[4, 5] = -length * sync_part.beta() ** 2 + rho * sx
        matrix[4, 5] *= get_dp_p_coeff(sync_part)
        return matrix

    def tilt(self, angle: float) -> np.ndarray:
        matrix = np.identity(7)
        matrix[0, 0] = matrix[1, 1] = +math.cos(angle)
        matrix[0, 2] = matrix[1, 3] = +math.sin(angle)
        matrix[2, 0] = matrix[3, 1] = -math.sin(angle)
        matrix[2, 2] = matrix[3, 3] = +math.cos(angle)
        return matrix

    def translation(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> np.ndarray:
        matrix = np.identity(7)
        matrix[0, 6] = x
        matrix[2, 6] = y
        matrix[4, 6] = z
        return matrix

    def kick(self, kx: float, ky: float, dE: float) -> np.ndarray:
        matrix = np.identity(7)
        matrix[1, 6] = kx
        matrix[3, 6] = ky
        matrix[5, 6] = dE
        return matrix

    def __call__(self, node: AccNode, sync_part: SyncParticle, part_index: int = 0) -> np.ndarray:
        if type(node) is DriftTEAPOT:
            length = node.getLength(part_index)
            return self.drift(length=length, sync_part=sync_part)
        elif type(node) is QuadTEAPOT:
            length = node.getLength(part_index)

            scale = 1.0
            if node.waveform:
                scale = node.waveform.getStrength()

            kq = scale * node.getParam("kq")
            return self.quad(length=length, kq=kq, sync_part=sync_part)

        elif type(node) is BendTEAPOT:
            nparts = node.getnParts()
            length = node.getLength(part_index)
            theta = node.getParam("theta") / (nparts - 1)
            return self.bend(length=length, theta=theta, sync_part=sync_part)

        elif type(node) is KickTEAPOT:
            nparts = node.getnParts()

            scale = 1.0
            if node.waveform:
                scale = node.waveform.getStrength()

            kx = scale * node.getParam("kx") / (nparts - 1)
            ky = scale * node.getParam("ky") / (nparts - 1)
            dE = node.getParam("dE") / (nparts - 1)
            return self.kick(kx, ky, dE)

        elif type(node) is TiltTEAPOT:
            angle = node.getTiltAngle()
            return self.tilt(angle)

        elif type(node) in self.ignore_node_types:
            return np.identity(7)

        else:
            if self.ignore_unknown:
                return self.drift(length=node.getLength(), sync_part=sync_part)
            raise NotImplementedError("Unsupported node type: {}. Set `ignore_unknown=True` to replace with drift.".format(type(node)))
