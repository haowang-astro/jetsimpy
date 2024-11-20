import numpy as np

# equal spacing
def Uniform(npoints):
    cells = np.linspace(0.0, np.pi, npoints)
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming forward-jet
def ForwardJetRes(theta_c, npoints):
    arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c), npoints)
    cells = np.sinh(arcsinhcells) * theta_c
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming counter-jet
def CounterJetRes(theta_c, npoints):
    cells = ForwardJetRes(theta_c, npoints)
    cells = np.pi - cells
    cells = np.flip(cells)
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming forward-jet & counter-jet
def ForwardCounterJetRes(theta_c, npoints):
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

# ------------- same to the above but the naming is not clear enough ------------ #

# cell edges assuming forward-jet
def NorthPole(theta_c, npoints):
    arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c), npoints)
    cells = np.sinh(arcsinhcells) * theta_c
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming counter-jet
def SouthPole(theta_c, npoints):
    cells = NorthPole(theta_c, npoints)
    cells = np.pi - cells
    cells = np.flip(cells)
    cells[0] = 0.0
    cells[-1] = np.pi
    return cells

# cell edges assuming forward-jet & counter-jet
def BothPoles(theta_c, npoints):
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