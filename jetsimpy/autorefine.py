import numpy as np
from . import jetsimpy_extension

class AutoGrid:
    def __init__(self, ntheta=128, jettype="forward"):
        self.ntheta = ntheta
        self.jettype = jettype
    
    def find_angle1(self, x, y):
        int1 = np.trapz(y * x ** 2 * np.sin(x), x)
        int2 = np.trapz(y * np.sin(x), x)
        return np.sqrt(2 * int1 / int2) if int2 != 0 else np.pi / 2
    
    def find_angle2(self, x, y1, y2):
        int1 = np.trapz(y1 * x ** 2 * np.sin(x), x)
        int2 = np.trapz(y1 * np.sin(x), x)
        th_c1 = np.sqrt(2 * int1 / int2) if int2 != 0 else np.pi / 2

        int1 = np.trapz(y2 * x ** 2 * np.sin(x), x)
        int2 = np.trapz(y2 * np.sin(x), x)
        th_c2 = np.sqrt(2 * int1 / int2) if int2 != 0 else np.pi / 2

        return min(th_c1, th_c2)
    
    def generate_cells_foward(self, xc):
        # construct grids
        arcsinhcells = np.linspace(0, np.arcsinh(np.pi / xc), self.ntheta)
        cells = np.sinh(arcsinhcells) * xc
        cells[-1] = np.pi

        return cells
    
    def generate_cells_both(self, xc1, xc2):
        # forward cells
        arcsinhcells = np.linspace(0, np.arcsinh(np.pi / 2 / xc1), int(self.ntheta / 2))
        cells1 = np.sinh(arcsinhcells) * xc1
        cells1[-1] = np.pi / 2

        # backward cells
        arcsinhcells = np.linspace(0, np.arcsinh(np.pi / 2 / xc2), int(self.ntheta / 2))
        cells2 = np.sinh(arcsinhcells) * xc2
        cells2[-1] = np.pi / 2
        cells2 = np.pi - np.flip(cells2)

        cells = np.hstack([cells1, cells2[1:]])
        return cells
    
    def grid1(self, x, y):
        if self.jettype == "forward":
            xc = self.find_angle1(x, y)
            cells = self.generate_cells_foward(xc)
        elif self.jettype == "both":
            xc1 = self.find_angle1(x[x < np.pi / 2], y[x < np.pi / 2])
            xc2 = self.find_angle1(np.flip(np.pi - x[x > np.pi / 2]), np.flip(y[x > np.pi / 2]))
            cells = self.generate_cells_both(xc1, xc2)
        else:
            raise RuntimeError("jettype can only be 'forward' or 'both'.")
        return cells

    def grid2(self, x, y1, y2):
        if self.jettype == "forward":
            xc = self.find_angle2(x, y1, y2)
            cells = self.generate_cells_foward(xc)
        elif self.jettype == "both":
            xc1 = self.find_angle2(x[x < np.pi / 2], y1[x < np.pi / 2], y2[x < np.pi / 2])
            xc2 = self.find_angle2(np.flip(np.pi - x[x > np.pi / 2]), np.flip(y1[x > np.pi / 2]), np.flip(y2[x > np.pi / 2]))
            cells = self.generate_cells_both(xc1, xc2)
        else:
            raise RuntimeError("jettype can only be 'forward' or 'both'.")
        return cells