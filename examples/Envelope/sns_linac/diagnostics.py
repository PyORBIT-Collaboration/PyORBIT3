import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis


class BunchMonitor:
    def __init__(self) -> None:
        self.twiss_calc = BunchTwissAnalysis()
        self.position_start = 0.0

        self.history = {}
        self.history["position"] = []
        self.history["rms_x"] = []
        self.history["rms_y"] = []
        self.history["rms_z"] = []
        self.history["kin_energy"] = []

    def __call__(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        node = params_dict["node"]
        position = params_dict["path_length"]

        if params_dict["old_pos"] == position:
            return
        if params_dict["old_pos"] + params_dict["pos_step"] > position:
            return
        params_dict["old_pos"] = position
        params_dict["count"] += 1

        sync_part = bunch.getSyncParticle()

        self.twiss_calc.analyzeBunch(bunch)

        cov_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                cov_matrix[i, j] = cov_matrix[j, i] = self.twiss_calc.getCorrelation(i, j)

        xrms = 1000.0 * np.sqrt(cov_matrix[0, 0])
        yrms = 1000.0 * np.sqrt(cov_matrix[2, 2])
        zrms = 1000.0 * np.sqrt(cov_matrix[4, 4])

        message = ""
        message += " s={:0.3f}".format(position + self.position_start)
        message += " xrms={:0.3f}".format(xrms)
        message += " yrms={:0.3f}".format(yrms)
        message += " zrms={:0.3f}".format(zrms)
        message += " node={}".format(node.getName())
        print(message)

        self.history["position"].append(position + self.position_start)
        self.history["rms_x"].append(xrms)
        self.history["rms_y"].append(yrms)
        self.history["rms_z"].append(zrms)
        self.history["kin_energy"].append(sync_part.kinEnergy())


