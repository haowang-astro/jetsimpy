import numpy as np
from numpy.typing import NDArray

# cell edges assuming forward-jet
def NorthPole(theta_c: float, npoints: int) -> NDArray:
    arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c), npoints)
    cells = np.sinh(arcsinhcells) * theta_c
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming counter-jet
def SouthPole(theta_c: float, npoints: int) -> NDArray:
    cells = NorthPole(theta_c, npoints)
    cells = np.pi - cells
    cells = np.flip(cells)
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming forward-jet & counter-jet
def BothPoles(theta_c: float, npoints: int) -> NDArray:
    half_points = int(npoints / 2) + 1
    arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c / 2.0), half_points)
    cells_n = np.sinh(arcsinhcells) * theta_c * 2.0
    if npoints % 2 == 0: # even numbers
        cells_n = cells_n / (2.0 - (cells_n[-1] - cells_n[-2]) / np.pi)
        cells_s = np.flip(np.pi - cells_n)
        cells = np.hstack([cells_n[:-1], cells_s[1:]])
    else:
        cells_n /= 2
        cells_s = np.flip(np.pi - cells_n)
        cells = np.hstack([cells_n, cells_s[1:]])
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# equal spacing
def Uniform(npoints: int) -> NDArray:
    cells = np.linspace(0.0, np.pi, npoints)
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells
