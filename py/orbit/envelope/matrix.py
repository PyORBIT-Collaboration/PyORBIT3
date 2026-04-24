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


class MatrixFactory:
    """Factory for 7x7 transfer matrices."""

    def __init__(self) -> None:
        self.ignore_node_types = [
            ApertureTEAPOT,
            BunchWrapTEAPOT,
            FringeFieldTEAPOT,
            MonitorTEAPOT,
            TurnCounterTEAPOT,
        ]

    @staticmethod
    def drift(length: float, gamma: float) -> np.ndarray:
        matrix = np.identity(7)
        matrix[0, 1] = length
        matrix[2, 3] = length
        matrix[4, 5] = length / gamma**2
        return matrix

    @staticmethod
    def quad(length: float, kq: float) -> np.ndarray:
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
        return matrix

    @staticmethod
    def bend(length: float, theta: float, gamma: float) -> np.ndarray:
        rho = length / theta

        cx = math.cos(theta)
        sx = math.sin(theta)

        matrix = np.identity(7)
        matrix[0, 0] = cx
        matrix[0, 1] = sx * rho
        matrix[0, 5] = (1.0 - cx) * rho
        matrix[1, 0] = -sx / rho
        matrix[1, 1] = cx
        matrix[1, 5] = sx
        matrix[2, 3] = length
        matrix[4, 0] = -sx
        matrix[4, 1] = -(1.0 - cx) * rho
        matrix[4, 5] = (length / gamma**2) - rho * (theta - sx)
        return matrix

    @staticmethod
    def tilt(angle: float) -> np.ndarray:
        cs = math.cos(angle)
        sn = math.sin(angle)

        matrix = np.identity(7)
        matrix[0, 0] = matrix[1, 1] = +cs
        matrix[0, 2] = matrix[1, 3] = +sn
        matrix[2, 0] = matrix[3, 1] = -sn
        matrix[2, 2] = matrix[3, 3] = +cs
        return matrix

    @staticmethod
    def kick(kx: float, ky: float, dE: float) -> np.ndarray:
        matrix = np.identity(7)
        matrix[1, 6] = kx
        matrix[3, 6] = ky
        matrix[5, 6] = dE
        return matrix
    
    @staticmethod
    def space_charge_2d(length: float, cov_matrix: np.ndarray, perveance: float) -> np.ndarray:
        # Start by assuming upright beam
        cx = 2.0 * math.sqrt(cov_matrix[0, 0])
        cy = 2.0 * math.sqrt(cov_matrix[2, 2])

        kappa_x = 2.0 * perveance / (cx * (cx + cy))
        kappa_y = 2.0 * perveance / (cy * (cx + cy))

        matrix = np.identity(7)
        matrix[1, 0] = kappa_x * length
        matrix[3, 2] = kappa_y * length
        return matrix
    
    @staticmethod
    def space_charge_3d(length: float, cov_matrix: np.ndarray, intensity: float) -> np.ndarray:
        raise NotImplementedError()

    def __call__(
        self, node: AccNode, sync_part: SyncParticle, part_index: int = 0
    ) -> np.ndarray:
        if type(node) is DriftTEAPOT:
            length = node.getLength(part_index)
            gamma = sync_part.gamma()
            return self.drift(length=length, gamma=gamma)

        elif type(node) is QuadTEAPOT:
            nparts = node.getnParts()
            length = node.getLength(part_index)

            scale = 1.0
            if node.waveform:
                scale = node.waveform.getStrength()

            kq = scale * node.getParam("kq")
            return self.quad(length=length, kq=kq)

        elif type(node) is BendTEAPOT:
            nparts = node.getnParts()
            length = node.getLength(part_index)
            theta = node.getParam("theta") / (nparts - 1)
            gamma = sync_part.gamma()
            return self.bend(length=length, theta=theta, gamma=gamma)

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
            raise NotImplementedError("Unsupported node type: {}".format(type(node)))
        