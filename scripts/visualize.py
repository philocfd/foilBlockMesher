import matplotlib.pyplot as plt
def visual_3D_point(x_coord, y_coord, z_coord):
    ax = plt.axes(projection="3d")
    ax.scatter3D(x_coord, y_coord, z_coord, c=z_coord, cmap="Greens")
    ax.set_xlabel("x axis")
    ax.set_ylabel("y axis")
    ax.set_zlabel("z axis")
    return ax