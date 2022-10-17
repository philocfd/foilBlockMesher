#!/usr/bin/env python
from pathlib import Path
from numpy import linspace, zeros, ones, sin, cos, arctan, pi
import matplotlib.pyplot as plt
import numpy as np
import argparse
def gen_blockmeshdict(foil_data, out_dir):
    """
    Write a `blockMeshDict` for a NACA foil at specified angle of attack.
    """
    # Foil geometry
    c = 1.0              # Geometric chord length
    # alpha = np.deg2rad(alpha_deg)  # Angle of attack (in radians)
    # NACA = [int(d) for d in foil]  # NACA 4-digit designation

    # Mesh dimensions
    scale = 1            # Scaling factor
    H = 8                # *Half* height of channel
    W = 0.5              # *Half* depth of foil (y-direction)
    D = 16               # Length of downstream section

    # Mesh resolution parameters
    Ni = 199             # Number of interpolation points along the foil
    Nx = 200          # Number of mesh cells along the foil
    ND = 150            # Number of cells in the downstream direction
    NT = 120         # Number of cells the transverse direction
    # Number of cells in the y-direction (along the foil axis)
    NW = 1

    # Expansion rates
    ExpT = 800           # Expansion rate in transverse direction
    ExpD = 80           # Expansion rate in the downstream direction
    ExpArc = 50          # Expansion rate along the inlet arc
    beta = linspace(0, pi, Ni)
    x = c*(0.5*(1 - cos(beta)))

    Xu = foil_data.upper_r[0, :].transpose()
    Zu = foil_data.upper_r[1, :].transpose()
    Xl = foil_data.lower_r[0, :].transpose()
    Zl = foil_data.lower_r[1, :].transpose()

    C_max_idx = np.where(foil_data.thickness == max(foil_data.thickness))[0][0]

    NoseX = (-H + Xu[C_max_idx])*cos(foil_data.aoa)
    NoseZ = -(-H + Xu[C_max_idx])*sin(foil_data.aoa)

    # Calculate the location of the vertices on the positive y-axis and put them in a matrix
    vertices = zeros((12, 3))

    vertices[0, :] = [NoseX[0], W, NoseZ[0]]
    vertices[1, :] = [Xu[C_max_idx], W, H]
    vertices[2, :] = [Xu[-1], W, H]
    vertices[3, :] = [D, W, H]
    vertices[4, :] = [Xu[0], W, Zu[0]]
    vertices[5, :] = [Xu[C_max_idx], W, Zu[C_max_idx]]
    vertices[6, :] = [Xl[C_max_idx], W, Zl[C_max_idx]]
    vertices[7, :] = [Xu[-1], W, Zu[-1]]
    vertices[8, :] = [D, W, Zu[-1]]
    vertices[9, :] = [Xl[C_max_idx], W, -H]
    vertices[10, :] = [Xu[-1], W, -H]
    vertices[11, :] = [D, W, -H]
    # annotations = ["P0", "P1", "P2", "P3", "P4", "P5",
    #                "P6", "P7", "P8", "P9", "P10", "P11"]
    # fig, ax = plt.subplots(figsize=(16, 9))
    # ax.scatter(vertices[:, 0], vertices[:, 2])
    # ax.set_xlabel("x axis")
    # ax.set_ylabel("z axis")
    # for i, label in enumerate(annotations):
    #     plt.annotate(label, (vertices[:, 0][i], vertices[:, 2][i]))
    # plt.show()
    # Create vertices for other side (negative y-axis)
    vertices2 = vertices.copy()
    vertices2[:, 1] *= -1
    vertices = np.vstack((vertices, vertices2))

    # Edge 4-5 and 16-17
    pts1 = np.concatenate([Xu[1:C_max_idx], W*ones(np.shape(Xu[1:C_max_idx])),
                           Zu[1:C_max_idx]], axis=1)
    pts5 = np.concatenate([pts1[:, 0], -pts1[:, 1], pts1[:, 2]], axis=1)

    # Edge 5-7 and 17-19
    pts2 = np.concatenate([Xu[C_max_idx + 1:Ni - 1],
                           W*ones(np.shape(Xu[C_max_idx + 1:Ni - 1])),
                           Zu[C_max_idx + 1:Ni - 1]], axis=1)
    pts6 = np.concatenate([pts2[:, 0], -pts2[:, 1], pts2[:, 2]], axis=1)

    # Edge 4-6 and 16-18
    pts3 = np.concatenate([Xl[1:C_max_idx], W*ones(np.shape(Xl[1:C_max_idx])),
                           Zl[1:C_max_idx]], axis=1)
    pts7 = np.concatenate([pts3[:, 0], -pts3[:, 1], pts3[:, 2]], axis=1)

    # Edge 6-7 and 18-19
    pts4 = np.concatenate([Xl[C_max_idx + 1:Ni - 1],
                          W*ones(np.shape(Xl[C_max_idx + 1:Ni - 1])),
                          Zl[C_max_idx + 1:Ni - 1]], axis=1)
    pts8 = np.concatenate([pts4[:, 0], -pts4[:, 1], pts4[:, 2]], axis=1)

    # Edge 0-1 and 12-13
    pts9 = np.array([-H*cos(pi/4) + Xu[C_max_idx, 0], W, H*sin(pi/4)])
    pts11 = np.array([pts9[0], -pts9[1], pts9[2]])

    # Edge 0-9 and 12-21
    pts10 = np.array([-H*cos(pi/4) + Xu[C_max_idx, 0], W, -H*sin(pi/4)])
    pts12 = np.array([pts10[0], -pts10[1], pts10[2]])

    # Calculate number of mesh points along 4-5 and 4-6
    Nleading = int((x[C_max_idx]/c)*Nx)

    # Calculate number of mesh points along 5-7 and 6-7
    Ntrailing = Nx - Nleading

    # Open file
    out_dir = Path(out_dir)
    f = open(out_dir/"blockMeshDict", "w")

    # Write file
    f.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n")
    f.write("| =========                 |                                                 | \n")
    f.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
    f.write("|  \\\\    /   O peration     | Version:  3.0.x                                 | \n")
    f.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      | \n")
    f.write("|    \\\\/     M anipulation  |                                                 | \n")
    f.write("\\*---------------------------------------------------------------------------*/ \n")
    f.write("FoamFile                                                                        \n")
    f.write(
        "{                                                                               \n")
    f.write("    version     2.0;                                                            \n")
    f.write("    format      ascii;                                                          \n")
    f.write("    class       dictionary;                                                     \n")
    f.write("    object      blockMeshDict;                                                  \n")
    f.write("}                                                                               \n")
    f.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n")
    f.write("\n")
    f.write("convertToMeters %f; \n" % scale)
    f.write("\n")
    f.write("vertices \n")
    f.write("( \n")
    for vertex in vertices:
        f.write("    (%f %f %f)\n" % tuple(vertex))
    f.write("); \n")
    f.write("\n")
    f.write("blocks \n")
    f.write("( \n")
    f.write("    hex (4 5 1 0 16 17 13 12)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n" % (
        Nleading, NT, NW, 1/ExpArc, 1/ExpArc, ExpT, ExpT, ExpT, ExpT))
    f.write("    hex (5 7 2 1 17 19 14 13)     (%i %i %i) simpleGrading (1 %f 1) \n" % (
        Ntrailing, NT, NW, ExpT))
    f.write("    hex (7 8 3 2 19 20 15 14)     (%i %i %i) simpleGrading (%f %f 1) \n" % (
        ND, NT, NW, ExpD, ExpT))
    f.write("    hex (16 18 21 12 4 6 9 0)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n" % (
        Nleading, NT, NW, 1/ExpArc, 1/ExpArc, ExpT, ExpT, ExpT, ExpT))
    f.write("    hex (18 19 22 21 6 7 10 9)    (%i %i %i) simpleGrading (1 %f 1) \n" % (
        Ntrailing, NT, NW, ExpT))
    f.write("    hex (19 20 23 22 7 8 11 10)   (%i %i %i) simpleGrading (%f %f 1) \n" % (
        ND, NT, NW, ExpD, ExpT))

    f.write("); \n")
    f.write("\n")
    f.write("edges \n")
    f.write("( \n")

    f.write("    spline 4 5 \n")
    f.write("        ( \n")
    for pt in np.array(pts1):
        f.write("            (%f %f %f) \n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 5 7 \n")
    f.write("        ( \n")
    for pt in np.array(pts2):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 4 6 \n")
    f.write("        ( \n")
    for pt in np.array(pts3):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 6 7 \n")
    f.write("        ( \n")
    for pt in np.array(pts4):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 16 17 \n")
    f.write("        ( \n")
    for pt in np.array(pts5):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 17 19 \n")
    f.write("        ( \n")
    for pt in np.array(pts6):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 16 18 \n")
    f.write("        ( \n")
    for pt in np.array(pts7):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 18 19 \n")
    f.write("        ( \n")
    for pt in np.array(pts8):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    arc 0 1 (%f %f %f) \n" % tuple(pts9))
    f.write("    arc 0 9 (%f %f %f) \n" % tuple(pts10))
    f.write("    arc 12 13 (%f %f %f) \n" % tuple(pts11))
    f.write("    arc 12 21 (%f %f %f) \n" % tuple(pts12))

    f.write("); \n")
    f.write("\n")
    f.write("boundary \n")
    f.write("( \n")

    f.write("    inlet \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (1 0 12 13) \n")
    f.write("            (0 9 21 12) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    outlet \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (11 8 20 23) \n")
    f.write("            (8 3 15 20) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    topAndBottom \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (3 2 14 15) \n")
    f.write("            (2 1 13 14) \n")
    f.write("            (9 10 22 21) \n")
    f.write("            (10 11 23 22) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    airfoil \n")
    f.write("    { \n")
    f.write("        type wall; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (5 4 16 17) \n")
    f.write("            (7 5 17 19) \n")
    f.write("            (4 6 18 16) \n")
    f.write("            (6 7 19 18) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("); \n")
    f.write(" \n")
    f.write("mergePatchPairs \n")
    f.write("( \n")
    f.write("); \n")
    f.write(" \n")
    f.write("// ************************************************************************* // \n")

    # Close file
    f.close()