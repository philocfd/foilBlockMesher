from scripts.blockmeshdict import dufoil, gen_blockmeshdict
import pandas as pd
import matplotlib.pyplot as plt
foil = dufoil(
    r"C:\Users\33937\Dropbox\Code\CAD\curiosityFluidsAirfoilMesher\example\DU30_A17_coords.dat", aoa=20)
print(foil.circum, foil.area)
# upper_aft_r, lower_aft_r = foil.rotate_foil()
# xu = upper_aft_r[0, :].conj().transpose()
# print(xu, type(xu), xu.shape)
# xu1 = upper_aft_r[0, :].transpose()
# print(xu1, type(xu1), xu1.shape)
# assert xu.all() == xu1.all()
# upper_aft_r_df = pd.DataFrame(upper_aft_r)
# # print(upper_aft_r_df.head())
# foil.visualize(upper_aft_r_df.transpose())
# foil.visualize_upper()
# fig, ax = plt.subplots()
# ax.plot(foil.upper.iloc[:, 0], foil.thickness)
# plt.show()
# gen_blockmeshdict(foil_data=foil)
foil.write_coord(
    r"D:\Work\2022\06June\wall turbulence\blockMesh\DU30Re1M\AOA_20")
