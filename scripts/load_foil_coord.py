import pandas as pd
import matplotlib.pyplot as plt
foil_coord = pd.read_table(r"../foil_coord/DU30_A17_coords.txt", sep="\s+", skiprows=8, header=None)
fig, ax = plt.subplots(figsize=(16, 9))
ax.scatter(foil_coord.iloc[:,0].values, foil_coord.iloc[:, 1].values)
ax.set_xlabel("x axis")
ax.set_ylabel("z axis")
plt.show()
