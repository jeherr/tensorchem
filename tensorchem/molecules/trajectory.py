"""
Trajectories hold Geometries from a method which generates multiple
Geometries which may be connected through sampling dynamics, an optimization
trajectory, or by applying perturbations (e.g. normal mode sampling)
"""

import os

from tensorchem.molecules.geometry import Geometry


class Trajectory:
    def __init__(self):
        self.geometries = []
        return
    
    def to_json(self):
        data_dict = {"geometries": [geom.to_json() for geom in self.geometries]}
        return data_dict

    @classmethod
    def from_json(cls, json_data: dict):
        new_traj = cls()
        new_traj.geometries = [Geometry.from_json(geom) for geom in json_data["geometries"]]
        return new_traj

    def write_xyz_trajectory(self, filename: str, filepath: str = None):
        if filepath is None:
            filepath = os.getcwd()
        with open(os.path.join(filepath, filename), "w") as f:
            for geom in self.geometries:
                f.write(str(geom))
                f.write("\n")


class OptTrajectory(Trajectory):
    def __init__(self):
        super().__init__()
        self.opt_algo = None
        return

    def to_json(self):
        data_dict = super(OptTrajectory, self).to_json()
        data_dict["opt_algorithm"] = self.opt_algo
        return data_dict

    def from_json(self, json_data: dict):
        new_geom = super(OptTrajectory, self).from_json(json_data)
        new_geom.opt_algo = json_data["opt_algo"]
        return new_geom


class NMSTrajectory(Trajectory):
    def __init__(self):
        super().__init__()
        self.nms_temp = None
        self.conformer = None
        return

    def to_json(self):
        data_dict = super(NMSTrajectory, self).to_json()
        data_dict['nms_temp'] = self.nms_temp
        data_dict['conformer'] = self.conformer.to_json()
        return data_dict

    def from_json(self, json_data: dict):
        new_geom = super(NMSTrajectory, self).from_json(json_data)
        new_geom.nms_temp = json_data["nms_temp"]
        new_geom.conformer = Geometry.from_json(json_data["conformer"])
        return new_geom
