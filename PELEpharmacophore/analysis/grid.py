import glob
import numpy as np
import os
import re
import PELEpharmacophore.helpers as hl

class Grid:

    def __init__(self, center, radius):
        self.center = center  # user-defined for WaterMap, otherwise ligand COI
        self.radius = radius  # user-defined for WaterMap, otherwise max coordinates of the ligand
        self.side = int(radius) * 2
        self.n_voxels = int(self.side) ** 3
        self.active_voxels = []

    def generate_voxels(self):
        center = self.center
        side = self.side
        self.voxels = []

        print("Creating {} voxels.".format(self.n_voxels))

        self.v1 = np.add(center, [-side / 2, -side / 2, -side / 2])  # lowest coordinate vertex
        self.v8 = np.add(center, [side / 2, side / 2, side / 2])  # highest coordinate vertex

        first_voxel_center = self.v1 + np.array([0.5, 0.5, 0.5])
        last_voxel_center = self.v8 + np.array([-0.5, -0.5, -0.5])

        samples = self.side

        x = np.linspace(first_voxel_center[0], last_voxel_center[0], samples)
        y = np.linspace(first_voxel_center[1], last_voxel_center[1], samples)
        z = np.linspace(first_voxel_center[2], last_voxel_center[2], samples)

        x_coord, y_coord, z_coord = np.meshgrid(x, y, z)
        coordinate_grid = np.array([x_coord, y_coord, z_coord])

        voxels = list(zip(coordinate_grid[0, :, :].reshape(1, self.n_voxels).tolist()[0],
                     coordinate_grid[1, :, :].reshape(1, self.n_voxels).tolist()[0],
                     coordinate_grid[2, :, :].reshape(1, self.n_voxels).tolist()[0]))

        for v in voxels:
            voxel = Voxel(v)
            self.voxels.append(voxel)

    def add_active_voxel(self, voxel):
        self.active_voxels.append(voxel)


class Voxel:

    def __init__(self, v):
        self.center = v
        self.freq_dict = None
        self.origin_dict = None

    def count_feature(self, feature):
        self.freq_dict = hl.frequency_dict(self.freq_dict, feature, 1)

    def add_origin(self, feature, origin):
        self.origin_dict = hl.list_dict(self.origin_dict, feature, origin)
