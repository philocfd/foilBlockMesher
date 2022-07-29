import numpy as np


class foil():
    def __init__(self, foil="0012", alpha_deg=4):
        self.foil_type = foil
        # convert aoa in degree to radians
        self.aoa = np.deg2rad(alpha_deg)
        self.chordL = 1.0
        self.NACA = [int(d) for d in foil]
