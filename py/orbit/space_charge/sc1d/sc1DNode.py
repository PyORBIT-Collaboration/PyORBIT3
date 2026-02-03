"""
Module. Includes classes for 1D longidutinal space charge accelerator nodes.
"""

import sys
import os
import math

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import LSpaceChargeCalc
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNodeBunchTracker
from orbit.utils import consts
from orbit.utils import orbitFinalize
from orbit.teapot import DriftTEAPOT


class SC1D_AccNode(DriftTEAPOT):
    """Longitudinal space charge node."""

    def __init__(
        self,
        b_a: float,
        phase_length: float,
        nmacros_min: float,
        use_sc: float,
        nbins: float,
        nmodes: int = None,
        use_grad: bool = False,
        name="long sc node",
    ) -> None:
        """
        Constructor. Creates the SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(b_a, phase_length, nmacros_min, use_sc, nbins)
        self.setNumModes(nmodes)
        # self.setUseGrad(use_grad)
        self.setType("long sc node")
        self.setLength(0.0)

    def setUseGrad(self, use_grad: bool) -> None:
        """Sets whether to use gradient-based solver instead of impedance solver."""
        self.lspacecharge.setUseGrad(int(use_grad))

    def setNumModes(self, n: int) -> None:
        """Sets number of FFT modes used to calculate energy kick."""
        self.lspacecharge.setNumModes(n)

    def trackBunch(self, bunch: Bunch) -> None:
        """
        The SC1D-teapot class implementation of the
        AccNodeBunchTracker class trackBunch(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        self.lspacecharge.trackBunch(bunch)  # track method goes here

    def track(self, params_dict: dict) -> None:
        """
        The SC1D-teapot class implementation of the
        AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = params_dict["bunch"]
        self.lspacecharge.trackBunch(bunch)

    def assignImpedance(self, py_cmplx_arr: list[float]) -> None:
        self.lspacecharge.assignImpedance(py_cmplx_arr)


class FreqDep_SC1D_AccNode(DriftTEAPOT):
    """Longitudinal space charge node (frequency-dependent)."""

    def __init__(
        self,
        b_a: float,
        phase_length: float,
        nmacros_min: int,
        use_sc: int,
        nbins: int,
        bunch: Bunch,
        imp_dict: dict,
        name: str = "freq. dep. long sc node",
    ) -> None:
        """
        Constructor. Creates the FreqDep_SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(
            b_a, phase_length, nmacros_min, use_sc, nbins
        )
        self.setType("freq. dep. long sc node")
        self.setLength(0.0)
        self.phase_length = phase_length
        self.nbins = nbins
        self.localDict = imp_dict
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = len(self.freq_tuple) - 1
        self.z_tuple = self.localDict["z_imp"]
        self.c = consts.speed_of_light
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins // 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = interp(freq_mode, self.freq_range, self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)

    def trackBunch(self, bunch: Bunch) -> None:
        """
        The FreqDep_SC1D-teapot class implementation of
        the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins // 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = interp(freq_mode, self.freq_range, self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)

    def track(self, params_dict: dict) -> None:
        """
        The FreqDep_SC1D-teapot class implementation of
        the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = params_dict["bunch"]
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins // 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = interp(freq_mode, self.freq_range, self.freq_tuple, self.z_tuple)
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)


class BetFreqDep_SC1D_AccNode(DriftTEAPOT):
    """Longitudinal space charge node (frequency- and velocity-dependent)."""

    def __init__(
        self,
        b_a: float,
        phase_length: float,
        nmacros_min: float,
        use_sc: int,
        nbins: int,
        bunch: Bunch,
        imp_dict: dict,
        name: str = "freq. dep. long sc node",
    ) -> None:
        """
        Constructor. Creates the BetFreqDep_SC1D-teapot element.
        """
        DriftTEAPOT.__init__(self, name)
        self.lspacecharge = LSpaceChargeCalc(
            b_a, phase_length, nmacros_min, use_sc, nbins
        )
        self.setType("beta-freq. dep. long sc node")
        self.setLength(0.0)
        self.phase_length = phase_length
        self.nbins = nbins
        self.localDict = imp_dict
        self.bet_tuple = self.localDict["betas"]
        self.bet_range = len(self.bet_tuple) - 1
        self.freq_tuple = self.localDict["freqs"]
        self.freq_range = len(self.freq_tuple) - 1
        self.z_bf = self.localDict["z_imp"]
        self.c = consts.speed_of_light
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(
                BetaRel,
                freq_mode,
                self.bet_range,
                self.freq_range,
                self.bet_tuple,
                self.freq_tuple,
                self.z_bf,
            )
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)

    def trackBunch(self, bunch: Bunch) -> None:
        """
        The BetFreqDep_SC1D-teapot class implementation of
        the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(
                BetaRel,
                freq_mode,
                self.bet_range,
                self.freq_range,
                self.bet_tuple,
                self.freq_tuple,
                self.z_bf,
            )
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)

    def track(self, params_dict: dict) -> None:
        """
        The BetFreqDep_SC1D-teapot class implementation of
        the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = params_dict["bunch"]
        BetaRel = bunch.getSyncParticle().beta()
        Freq0 = (BetaRel * self.c) / self.phase_length
        Z = []
        for n in range(self.nbins / 2 - 1):
            freq_mode = Freq0 * (n + 1)
            z_mode = bilinterp(
                BetaRel,
                freq_mode,
                self.bet_range,
                self.freq_range,
                self.bet_tuple,
                self.freq_tuple,
                self.z_bf,
            )
            Z.append(z_mode)
        self.lspacecharge.assignImpedance(Z)
        self.lspacecharge.trackBunch(bunch)


def interp(x: float, n_tuple: int, x_tuple: list[float], y_tuple: list[float]) -> float:
    """
    Linear interpolation: Given n-tuple + 1 points,
    x_tuple and y_tuple, routine finds y = y_tuple
    at x in x_tuple. Assumes x_tuple is increasing array.
    """
    if x < x_tuple[0]:
        y = y_tuple[0]
        return y
    if x > x_tuple[n_tuple]:
        y = y_tuple[n_tuple]
        return y
    dxp = x - x_tuple[0]
    for n in range(n_tuple):
        dxm = dxp
        dxp = x - x_tuple[n + 1]
        dxmp = dxm * dxp
        if dxmp <= 0:
            break
    y = (-dxp * y_tuple[n] + dxm * y_tuple[n + 1]) / (dxm - dxp)
    return y


def bilinterp(
    x: float,
    y: float,
    nx_tuple: int,
    ny_tuple: int,
    x_tuple: list[float],
    y_tuple: list[float],
    fxy: list[list[float]],
) -> float:
    """
    Bilinear interpolation: Given nx-tuple + 1 x-points,
    ny-tuple + 1 y-points, x_tuple and y_tuple,
    routine finds f(x, y) = fxy at (x, y) in (x_tuple, y_tuple).
    Assumes x_tuple and y_tuple are increasing arrays.
    """
    f_tuple = []
    if x < x_tuple[0]:
        for ny in range(ny_tuple + 1):
            vf = fxy[0][ny]
            f_tuple.append(vf)
    elif x > x_tuple[nx_tuple]:
        for ny in range(ny_tuple + 1):
            vf = fxy[x_tuple][ny]
            f_tuple.append(vf)
    else:
        dxp = x - x_tuple[0]
        for nx in range(nx_tuple):
            dxm = dxp
            dxp = x - x_tuple[nx + 1]
            dxmp = dxm * dxp
            if dxmp <= 0:
                break
    for ny in range(ny_tuple + 1):
        vf = (-dxp * fxy[nx][ny] + dxm * fxy[nx + 1][ny]) / (dxm - dxp)
        f_tuple.append(vf)
    f = interp(y, ny_tuple, y_tuple, f_tuple)
    return f
