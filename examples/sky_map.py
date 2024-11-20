# sky map of a gaussian jet

import numpy as np
from matplotlib import pyplot as plt
import jetsimpy
from matplotlib.patches import Ellipse

# put the parameters in a dictionary
P = dict(
    Eiso = 1e52,        # (Jet) Isotropic equivalent energy
    lf = 300,           # (Jet) Lorentz factor
    theta_c = 0.1,      # (Jet) half opening angle
    n0 = 1,             # (ISM) constant number density
    A = 0,              # (ISM) wind amplitude
    eps_e = 0.1,        # (Radiation) epsilon_e
    eps_b = 0.01,       # (Radiation) epsilon_b
    p = 2.17,           # (Radiation) electron power index
    theta_v = 0.4,      # (Radiation) viewing angle
    d = 474.33,         # (radiation) distance (Mpc)
    z = 0.1,            # (radiation) redshift
)

# ---------- (step 1) hydro simulation of the jet ---------- #

# solve jet
jet = jetsimpy.Jet(
    jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    P["A"],          # scale of wind density
    P["n0"],         # constant number density
    grid=jetsimpy.ForwardJetRes(P["theta_c"], 129),   # resolution
)

# ---------- (step 2) calculate intensity at some pixels ---------- #

# time and frequency of sky map
tday = 10
nu = 3e9
tsecond = tday * 3600 * 24

# make an image with half width: (centroid + 3 * size_x) [unit: mas],
# such that the whole jet is visible in the image
offset = jet.Offset(tsecond, nu, P)
size_x = jet.SizeX(tsecond, nu, P)
size_y = jet.SizeY(tsecond, nu, P)    # why? just for fun!
half_width = offset + 3 * size_x

# intensity map data matrix: resolution 300 x 300
resolution = 300
x_tilde = np.linspace(- half_width, half_width, resolution)    # x pixel coordinates
y_tilde = np.linspace(- half_width, half_width, resolution)    # y pixel coordinates
X, Y = np.meshgrid(x_tilde, y_tilde)
sky_map = jet.IntensityOfPixel(tsecond, nu, X, Y, P)

# configurate the plot
fig, ax = plt.subplots()
Inorm = sky_map / sky_map.max()    # normalize to 1
Inorm = np.rot90(Inorm, 2)         # rotate to align the plt.imshow() extent
im = ax.imshow(Inorm, interpolation='gaussian', cmap="inferno", extent=[-half_width, half_width, -half_width, half_width])
im.set_clim(vmin=0.0, vmax=1.0)

# plot the origin, centroid, and [optional] size contour just for fun! (you are welcome to uncomment the line below.)
ax.scatter(0, 0, marker="*", color="white", label="origin")
ax.scatter(offset, 0.0, marker="+", color="white", label="centroid")
#ax.add_patch(Ellipse((offset, 0), size_x * 2, size_y * 2, edgecolor="white", fill=False, linestyle="--", linewidth=1))

# show the image!
ax.set_xlabel(r"$\tilde{x}$ [mas]", fontsize=12)
ax.set_ylabel(r"$\tilde{y}$ [mas]", fontsize=12)
ax.legend(loc="upper right", frameon=False, labelcolor="white", fontsize=12, handletextpad=0.1)
fig.tight_layout()
plt.show()