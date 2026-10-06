import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis


class BunchMonitor:
    def __init__(self) -> None:
        self.twiss_calc = BunchTwissAnalysis()
        self.position_start = 0.0

        keys = [
            "s",
            "kin_energy",
            "gamma",
            "beta",
            "rms_x",
            "rms_y",
            "rms_z",
            "eps_x",
            "eps_x_n",
            "eps_y",
            "eps_y_n",
            "eps_z",
            "eps_z_n",
            "eps_1",
            "eps_2",
            "eps_3",
        ]
        self.history = {key: [] for key in keys}

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
                cov_matrix[i, j] = self.twiss_calc.getCorrelation(i, j)
                cov_matrix[j, i] = cov_matrix[i, j]

        centroid = np.zeros(6)
        for i in range(6):
            centroid[i] = self.twiss_calc.getAverage(i)

        results = {}
        results["s"] = position
        results["gamma"] = sync_part.gamma()
        results["beta"] = sync_part.beta()
        results["kin_energy"] = sync_part.kinEnergy()
        results["rms_x"] = np.sqrt(cov_matrix[0, 0])
        results["rms_y"] = np.sqrt(cov_matrix[2, 2])
        results["rms_z"] = np.sqrt(cov_matrix[4, 4])
        results["eps_x"] = np.sqrt(np.linalg.det(cov_matrix[0:2, 0:2]))
        results["eps_y"] = np.sqrt(np.linalg.det(cov_matrix[2:4, 2:4]))
        results["eps_z"] = np.sqrt(np.linalg.det(cov_matrix[4:6, 4:6]))

        poisson_matrix = np.zeros_like(cov_matrix)
        for i in range(0, 6, 2):
            poisson_matrix[i, i + 1] = +1.0
            poisson_matrix[i + 1, i] = -1.0
        eigvals = np.linalg.eigvals(cov_matrix @ poisson_matrix)
        eigvals = np.imag(eigvals)
        eigvals = eigvals[eigvals > 0]
        results["eps_1"] = eigvals[0]
        results["eps_2"] = eigvals[1]
        results["eps_3"] = eigvals[2]

        results["eps_x_n"] = results["eps_x"] * results["gamma"] * results["beta"]
        results["eps_y_n"] = results["eps_y"] * results["gamma"] * results["beta"]
        results["eps_z_n"] = results["eps_z"] / results["beta"]

        for key in results:
            self.history[key].append(results[key])

        message = ""
        message += " s={:0.3f}".format(position + self.position_start)
        message += " rms_x={:0.3f}".format(1000.0 * results["rms_x"])
        message += " rms_y={:0.3f}".format(1000.0 * results["rms_y"])
        message += " rms_z={:0.3f}".format(1000.0 * results["rms_z"])
        message += " node={}".format(node.getName())
        print(message)
