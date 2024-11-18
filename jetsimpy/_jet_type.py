import numpy as np

# top-hat jet
def TopHat(theta_c, Eiso, lf0=1e100):
    theta = theta = np.linspace(0, np.pi, 10000)
    
    energy = Eiso * np.ones_like(theta)
    energy[theta > theta_c] = 0
    
    lf = np.ones_like(theta)
    lf[theta <= theta_c] = lf0

    return (theta, energy, lf)

# Gaussian jet
def Gaussian(theta_c, Eiso, lf0=1e100):
    theta = theta = np.linspace(0, np.pi, 10000)

    energy = Eiso * np.exp(- 0.5 * (theta / theta_c) ** 2)
    lf = (lf0 - 1) * np.exp(- 0.5 * (theta / theta_c) ** 2) + 1
    
    return (theta, energy, lf)

# power-law jet
def PowerLaw(theta_c, Eiso, lf0=1e100, s=4.0):
    theta = theta = np.linspace(0, np.pi, 10000)

    energy = Eiso * np.power(1 + (theta / theta_c) ** 2, - s / 2.0)
    lf = (lf0 - 1.0) * np.power(1 + (theta / theta_c) ** 2, - s / 2.0) + 1.0

    return (theta, energy, lf)