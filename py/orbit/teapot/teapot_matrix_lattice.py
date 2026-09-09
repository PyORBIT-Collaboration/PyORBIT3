"""
`TEAPOT_MATRIX_Lattice` is a subclass of `MAXTRIX_Lattice`. 

Matrices for each element are estimated from particle tracking data through
each node. It is assumed that RF cavities in the ring are of type
`RingRFTEAPOT`.
"""
import os
import math

from orbit.core.bunch import Bunch
from orbit.core.teapot_base import MatrixGenerator
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.matrix_lattice import BaseMATRIX
from orbit.matrix_lattice import MATRIX_Lattice
from orbit.teapot import BaseTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import RingRFTEAPOT
from orbit.utils import orbitFinalize


class TEAPOT_MATRIX_Lattice(MATRIX_Lattice):
    """Subclass of the `MATRIX_Lattice` class. Shell class for the BaseMATRIX nodes collection.
    
    A `TEAPOT_Lattice` and `Bunch` instance are needed to construct the transfer matrices
    for specified beam parameters.
    """
    def __init__(self, teapot_lattice: TEAPOT_Lattice, bunch: Bunch, name: str = None) -> None:
        MATRIX_Lattice.__init__(self, name)
        
        if not isinstance(teapot_lattice, TEAPOT_Lattice):
            orbitFinalize("Constructor orbit.teapot.TEAPOT_MATRIX_Lattice needs  the TEAPOT_Lattice instance.")
        
        if name is None:
            name = teapot_lattice.getName()
        
        self.setName(name)
        self.teapot_lattice = teapot_lattice
        
        self.bunch = Bunch()
        self.lost_bunch = Bunch()
        bunch.copyEmptyBunchTo(self.bunch)
        bunch.copyEmptyBunchTo(self.lost_bunch)
        
        self.matrix_generator = MatrixGenerator()

        def twissAction(params_dict: dict) -> None:
            """Build `MATRIX_Lattice` from `TEAPOT_Lattice`."""
            node = params_dict["node"]
            bunch = params_dict["bunch"]
            active_index = node.getActivePartIndex()
            n_parts = node.getnParts()
            length = node.getLength(active_index)
            if isinstance(node, BaseTEAPOT) and not isinstance(node, RingRFTEAPOT):
                self.matrix_generator.initBunch(bunch)
                node.track(params_dict)

                matrix_node = BaseMATRIX(node.getName() + "_" + str(active_index))
                matrix_node.addParam("matrix_parent_node", node)
                matrix_node.addParam("matrix_parent_node_type", node.getType())
                matrix_node.addParam("matrix_parent_node_n_nodes", n_parts)
                matrix_node.addParam("matrix_parent_node_active_index", active_index)
                matrix_node.setLength(length)
                
                self.matrix_generator.calculateMatrix(bunch, matrix_node.getMatrix())
                self.addNode(matrix_node)
            
            if isinstance(node, RingRFTEAPOT):
                rf_node = RingRFTEAPOT(node.getName())
                rf_node.setParamsDict(node.getParamsDict().copy())
                self.addNode(rf_node)

        action_container = AccActionsContainer()
        action_container.addAction(twissAction, AccActionsContainer.BODY)
        
        params_dict = {}
        params_dict["bunch"] = self.bunch
        params_dict["lostbunch"] = self.lost_bunch
        params_dict["position"] = 0.0
        params_dict["useCharge"] = self.teapot_lattice.getUseRealCharge()
        
        self.teapot_lattice.trackActions(action_container, params_dict)
        self.makeOneTurnMatrix()
        self.initialize()

    def getKinEnergy(self) -> float:
        return self.bunch.getSyncParticle().kinEnergy()

    def rebuild(self, kin_energy=-1.0):
        if kin_energy > 0.0:
            self.bunch.getSyncParticle().kinEnergy(kin_energy)

        for matrix_node in self.getNodes():
            if isinstance(matrix_node, BaseMATRIX):
                node = matrix_node.getParam("matrix_parent_node")
                active_index = matrix_node.getParam("matrix_parent_node_active_index")
                n_parts = matrix_node.getParam("matrix_parent_node_n_nodes")

                if n_parts != node.getnParts():
                    msg = " orbit.teapot.TEAPOT_MATRIX_Lattice class" + os.linesep
                    msg = msg + "  rebuild(kin_energy = -1.0) method" + os.linesep
                    msg = msg + "  TEAPOT node=" + node.getName() + os.linesep
                    msg = msg + "  has been changed!" + os.linesep
                    msg = msg + "  Stop!" + os.linesep
                    orbitFinalize(msg)

                self.matrix_generator.initBunch(self.bunch)

                params_dict = {}
                params_dict["bunch"] = self.bunch
                params_dict["node"] = node

                node.setActivePartIndex(active_index)
                node.track(params_dict)
                self.matrix_generator.calculateMatrix(self.bunch, matrix_node.getMatrix())
        self.makeOneTurnMatrix()

    def getRingParametersDict(self) -> dict[str, float]:
        """Returns dictionary of ring parameters from one-turn matrix.

        Overloads `getRingParametersDict(p,m)` method from `MATRIX_Lattice`.
        """
        momentum = self.bunch.getSyncParticle().momentum()
        mass = self.bunch.getSyncParticle().mass()
        return MATRIX_Lattice.getRingParametersDict(self, momentum, mass)

    def getRingMatrix(self) -> BaseMATRIX:
        """Returns one-turn matrix."""
        return MATRIX_Lattice.getOneTurnMatrix(self)

    def getRingOrbit(self, z0: list[float]) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Returns tuple ([(s, x),...],[(s, y),...] )."""
        return self.trackOrbit(z0)

    def getRingTwissDataX(self) -> tuple[list[tuple[float, float]], list[tuple[float, float]], list[tuple[float, float]]]:
        """Returns tuple ([(s, tune_x), ...], [(s, alpha_x),...], [(s, beta_x),...]).

        Overloads the `getRingTwissDataX(p,m)` method from `MATRIX_Lattice`.
        """
        res_dict = self.getRingParametersDict()
        alpha_x = res_dict["alpha x"]
        beta_x = res_dict["beta x [m]"]
        return self.trackTwissData(alpha_x, beta_x, "x")

    def getRingTwissDataY(self) -> tuple[list[tuple[float, float]], list[tuple[float, float]], list[tuple[float, float]]]:
        """Returns tuple ([(s, tune_y), ...], [(s, alpha_y),...], [(s, beta_y),...]).

        Overloads the `getRingTwissDataY(p,m)` method from `MATRIX_Lattice`.
        """
        res_dict = self.getRingParametersDict()
        alpha_y = res_dict["alpha y"]
        beta_y = res_dict["beta y [m]"]
        return self.trackTwissData(alpha_y, beta_y, "y")

    def getRingDispersionDataX(self) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Returns tuple  ([(s, disp_x), ...], [(s, disp_prime_x), ...]).

        Overloads `getRingDispersionDataX(p,m)` method from `MATRIX_Lattice`.
        """
        res_dict = self.getRingParametersDict()
        disp = res_dict["dispersion x [m]"]
        disp_p = res_dict["dispersion prime x"]
        momentum = res_dict["momentum [GeV/c]"]
        mass = res_dict["mass [GeV]"]
        return self.trackDispersionData(momentum, mass, disp, disp_p, "x")

    def getRingDispersionDataY(self) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Returns tuple  ([(s, disp_y), ...], [(s, disp_prime_y), ...]).

        Overloads `getRingDispersionDataX(p,m)` method from `MATRIX_Lattice`.
        """
        res_dict = self.getRingParametersDict()
        disp = res_dict["dispersion y [m]"]
        disp_p = res_dict["dispersion prime y"]
        momentum = res_dict["momentum [GeV/c]"]
        mass = res_dict["mass [GeV]"]
        return self.trackDispersionData(momentum, mass, disp, disp_p, "y")

    def getTransferTwissDataX(self, alpha_x: float, beta_x: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]], list[tuple[float, float]]]:
        """Propagates Twiss parameters (alpha_x, beta) through the lattice."""
        return self.trackTwissData(alpha_x, beta_x, "x")

    def getTransferTwissDataY(self, alpha_y: float, beta_y: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]], list[tuple[float, float]]]:
        """Propagates Twiss parameters (alpha_y, beta_y) through lattice."""
        return self.trackTwissData(alpha_y, beta_y, "y")

    def getTransferDispersionDataX(self, disp: float, disp_p: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Propagates x dispersion and dispersion prime (x) through lattice."""
        res_dict = self.getRingParametersDict()
        momentum = res_dict["momentum [GeV/c]"]
        mass = res_dict["mass [GeV]"]
        return self.trackDispersionData(momentum, mass, disp, disp_p, "x")

    def getTransferDispersionDataY(self, disp: float, disp_p: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Propagates y dispersion and dispersion prime through lattice."""
        res_dict = self.getRingParametersDict()
        momentum = res_dict["momentum [GeV/c]"]
        mass = res_dict["mass [GeV]"]
        return self.trackDispersionData(momentum, mass, disp, disp_p, "y")

    def getChromaticitiesXY(self) -> tuple[float, float]:
        """Calculates chromaticities for x and y planes."""
        (tune_x, tmp0, tmp1) = self.getRingTwissDataX()
        (tune_y, tmp0, tmp1) = self.getRingTwissDataY()
        tune_x = tune_x[1][-1]
        tune_y = tune_y[1][-1]
        self.matrix_generator.initBunchChromCoeff(self.bunch)

        def twissAction(params_dict: dict) -> None:
            node = params_dict["node"]
            if isinstance(node, BaseTEAPOT) == True and isinstance(node, RingRFTEAPOT) == False:
                node.track(params_dict)

        action_container = AccActionsContainer()
        action_container.addAction(twissAction, AccActionsContainer.BODY)

        params_dict = {}
        params_dict["bunch"] = self.bunch

        self.teapot_lattice.trackActions(action_container, params_dict)

        res_coeff = self.matrix_generator.calcChromCoeff(self.bunch)
        (coeff_x_dE, coeff_xp_dE, coeff_y_dE, coeff_yp_dE) = res_coeff

        momentum = self.bunch.getSyncParticle().momentum()
        mass = self.bunch.getSyncParticle().mass()
        kin_energy = self.bunch.getSyncParticle().kinEnergy()

        chrom_x = -(momentum / (2 * math.sin(2 * math.pi * tune_x))) * (coeff_x_dE + coeff_xp_dE)
        chrom_x = (momentum / (mass + kin_energy)) * chrom_x

        chrom_y = -(momentum / (2 * math.sin(2 * math.pi * tune_y))) * (coeff_y_dE + coeff_yp_dE)
        chrom_y = (momentum / (mass + kin_energy)) * chrom_y

        chrom_x = chrom_x / (2.0 * math.pi)
        chrom_y = chrom_y / (2.0 * math.pi)
        return (chrom_x, chrom_y)

    def getTunes(self) -> tuple[float, float]:
        """Return total phase advance / 2pi (not fractional tunes)"""
        twiss_data_x = self.getRingTwissDataX()
        twiss_data_y = self.getRingTwissDataY()

        tune_data_x = twiss_data_x[0]
        tune_data_y = twiss_data_y[0]

        _, tune_x = tune_data_x[-1]
        _, tune_y = tune_data_y[-1]
        return (tune_x, tune_y)