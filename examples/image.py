import numpy as np
from matplotlib import pyplot as plt
import jetsimpy
from matplotlib.patches import Ellipse

para = dict(
        n0 = 1,        # (Jet required) ism number density
        k = 0.0,          # (Jet required) wind power index
        E0 = 1e52,        # (Jet customized) Energy / LF parameter
        LF0 = 1000,         # (Jet customized) Energy / LF parameter
        theta_c = 0.1,    # (Jet customized) Energy / LF parameter
        z = 0.01,         # (radiation required) redshift
        d = 44,           # (radiation required) distance
        theta_v = 0.4,    # (radiation required) viewing angle
        eps_e = 0.01,      # (radiation customized) epsilon_e
        eps_b = 0.001,     # (radiation customized) epsilon_b
        p = 2.2,          # (radiation customized) electron power index
        xi_N = 1.0,       # (radiation customized) total fraction of accelerated electrons
        b = 0.0,
)

# generate initial data
E0 = para["E0"]
lf0 = para["LF0"]
theta_c = para["theta_c"]
theta = np.linspace(0, np.pi, 1000)
Eiso = E0 * np.exp(- 0.5 * (theta / theta_c) ** 2)
lf = (lf0 - 1.0) * np.exp(- 0.5 * (theta / theta_c) ** 2) + 1

theta_c_min = para["theta_c"]
arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c_min), 200)
cells = np.sinh(arcsinhcells) * theta_c_min
cells[-1] = np.pi

# calculate jetsimpy
afterglow = jetsimpy.Afterglow(
    theta,           
    Eiso,           # (required) energy distribution
    lf,               # (required) LF distribution
    0.0,
    para["n0"],
    spread=True,    # (optional) lateral expansion effect
    tmin=10,
    tmax=3.2e10,         # (optional) PDE maximum time
    grid=cells,
    coast=True
)

# setup figure
fig, axes = plt.subplots(3, 1, figsize = (4, 8), layout="constrained")

ts = [0.03, 0.5, 120]

ticks = [
    [-0.001, 0, 0.001],
    [-0.01, 0, 0.01],
    [-0.5, 0, 0.5]
]

for k in [0, 1, 2]:
    # solve image
    t = ts[k] * 3600 * 24
    nu = 3e9
    image = afterglow.solve_image(t, nu, inum=100, pnum=20, para=para)

    # setup image scales & criticle points
    width = image.half_width
    centroid = image.offset
    yscale = image.yscale
    xscale = image.xscale

    # orientation
    orientation = 2

    # intensity and polarization data
    I = np.array(image.intensity_image).T

    # intensity scale
    Inorm = I / I.max()
    Inorm = np.rot90(Inorm, orientation) # rotate

    # plot intensity
    im = axes[k].imshow(Inorm, interpolation='gaussian', cmap="inferno", origin='lower', extent=[-width,     width, -width, width])
    im.set_clim(vmin=0.0, vmax=1)

    # origin
    axes[k].scatter(0, 0, marker="*", color="white", label="origin")

    # offset and scale
    xc = centroid
    yc = 0.0
    sigmax = xscale
    sigmay = yscale
    axes[k].scatter(xc, yc, marker="+", color="white", label="centroid")
    axes[k].add_patch(Ellipse((xc, yc), sigmax * 2, sigmay * 2, edgecolor="white", fill=False,   linestyle="--", linewidth=1))

    # plot time indicator
    text_t = axes[k].text(- width * 0.9, width * 0.8, "t = " + str(round(t / 3600 / 24, 2)) + " d",   color="white", fontsize=10)

    # setup x/y labels
    axes[k].set_xlabel(r"$\tilde{x}$ [mas]", fontsize=12)
    axes[k].set_ylabel(r"$\tilde{y}$ [mas]", fontsize=12)
    axes[k].tick_params(axis="x",direction="in", pad=-15, color="white", labelcolor="white")
    axes[k].tick_params(axis="y",direction="in", pad=-18, color="white", labelcolor="white", labelrotation=90)
    axes[k].set_xticks(ticks[k], ticks[k])
    axes[k].set_yticks(ticks[k], ticks[k])

    # legend
    axes[k].legend(loc="upper right", frameon=False, labelcolor="white", fontsize=10, handletextpad=0.1)

## setup colorbar
clb = fig.colorbar(im, ax=axes[0], location="top", fraction=0.05)
clb.set_label(r"$I_{\nu}\ /\ I_{\nu}^{\rm max}$", fontsize=15)
plt.show()