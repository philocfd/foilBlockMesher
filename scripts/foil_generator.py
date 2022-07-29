from __future__ import annotations, division, print_function
import numpy as np
import os
from CAD.foilBlockMesher.scripts.visualize import visual_3D_point
import matplotlib.pyplot as plt
import pandas as pd
import math
from pathlib import Path


def arc_length(points):
    """
    Calculate the lenght of a poly line made of a lot of points
    """
    points_array = points.to_numpy()
    length = 0
    list_front = points_array[:-1, :]
    list_back = points_array[1:, :]
    for point_front, point_back in zip(list_front, list_back):
        length = length+math.dist(point_front, point_back)
    return length


class nacafoil():
    def __init__(self, foil="0012", alpha_deg=4) -> None:
        self.foil_type = foil
        self.aoa = np.deg2rad(alpha_deg)  # Angle of attack on radians
        self.c = 1.0  # Chord length of foil
        self.NACA = [int(d) for d in foil]


class dufoil():
    """
    Class handling du foil geometry coordinates
    """

    def __init__(self, data_file, aoa=0):
        self.coord = pd.read_table(
            data_file, sep="\s+", skiprows=1, header=None)
        # To make the pressure center coincide with the coordinate origin
        self.coord.iloc[:, 0] = self.coord.iloc[:, 0]-0.25
        # dataFrame of upper and lower surface
        self.upper = self.coord.iloc[:int(len(self.coord.index)/2)]
        self.lower = self.coord.iloc[int(len(self.coord.index)/2):]
        self.chordl = 1.0
        self.aoa = np.deg2rad(aoa)
        # coordinate after rotation
        self.upper_r, self.lower_r = self.rotate_foil()
        # To make the direction pointing to the positive x axis
        self.upper_r = self.upper_r[:, ::-1]

        self.thickness = self.upper.iloc[:, 1].values - \
            self.lower.iloc[:, 1].values
        # circumference of foil
        self.circum = self.circum()
        # airea of foil
        self.area = np.sum(self.thickness)/self.thickness.size

    def rotate_foil(self):
        """
        Rotate the foil coordinate according to the value of angle of attack.
        Rotation center is (0,0)
        Rotation angle is positive when rotated counterclockwise

        return:
        2*n array of foil coordinate
        """
        R = np.matrix([[np.cos(self.aoa), np.sin(self.aoa)],
                       [-np.sin(self.aoa), np.cos(self.aoa)]])
        upper = R*np.vstack((self.upper.iloc[:, 0].values,
                             self.upper.iloc[:, 1].values))
        lower = R*np.vstack((self.lower.iloc[:, 0].values,
                             self.lower.iloc[:, 1].values))
        return upper, lower

    def visualize(self, data):
        """
        visualize foil coordinate data
        data: dataFrame for coordinate
        """
        fig, ax = plt.subplots(figsize=(9, 18))
        ax.scatter(data.iloc[:, 0].values,
                   data.iloc[:, 1].values)
        ax.set_xlabel("x axis")
        ax.set_ylabel("z axis")
        ax.set_xlim(min(data.iloc[:, 0].values)-0.25, 1)
        ax.set_ylim(-0.5, 0.5)
        plt.show()

    def visualize_upper(self):
        self.visualize(self.upper)

    def visualize_lower(self):
        self.visualize(self.lower)

    def circum(self):
        upper_length = arc_length(self.upper)
        lower_lenth = arc_length(self.lower)
        return upper_length+lower_lenth

    def write_coord(self, out_put_dir):
        "write coordinates for sample in OpenFOAM"
        with open(out_put_dir + os.sep+"upper_coord.dat", "w") as f:
            for row in self.upper_r.T:
                f.write("( {}\t{}\t{} )\n".format(row[0, 0], row[0, 1], 0))
        with open(out_put_dir + os.sep+"lower_coord.dat", "w") as f:
            for row in self.lower_r.T:
                f.write("( {}\t{}\t{} )\n".format(row[0, 0], row[0, 1], 0))
