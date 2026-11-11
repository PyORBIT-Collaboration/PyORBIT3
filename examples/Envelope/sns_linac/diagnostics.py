from orbit.envelope import Envelope
from orbit.envelope import EnvelopeAnalysis


class BunchMonitor:
    def __init__(self) -> None:
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

        analysis = EnvelopeAnalysis(Envelope(bunch=bunch))
        results = {"s": position, **analysis.to_dict()}
        for key in self.history:
            self.history[key].append(results[key])

        message = ""
        message += " s={:0.3f}".format(position + self.position_start)
        message += " rms_x={:0.3f}".format(1000.0 * results["rms_x"])
        message += " rms_y={:0.3f}".format(1000.0 * results["rms_y"])
        message += " rms_z={:0.3f}".format(1000.0 * results["rms_z"])
        message += " node={}".format(node.getName())
        print(message)
