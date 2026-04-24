import math

import numpy as np


class MatrixFactory:
    """Factory for 7x7 transfer matrices."""

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
        if (kq > 0):
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
        elif (kq < 0):
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

        cx  = math.cos(theta)
        sx  = math.sin(theta)
        
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
    
    