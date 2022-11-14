#!/usr/bin/env python
"""
multipath_wvfrm.py

Compute all eigenrays for a source-receiver pair
    and use the weakly non-linear waveform methods
    to compute each contribution to the waveform

Philip Blom (pblom@lanl.gov)

"""

import sys
import os

import numpy as np

import matplotlib.pyplot as plt 
import matplotlib.cm as cm

from mpl_toolkits.axes_grid1 import make_axes_locatable

def run(specification, spec_format='zTuvdp', alt_max=None):

    atmo = np.loadtxt(specification)
    z = atmo[:, spec_format.find('z')]
    u = atmo[:, spec_format.find('u')]
    v = atmo[:, spec_format.find('v')]
    c = np.sqrt(0.14 * atmo[:, spec_format.find('p')] / atmo[:, spec_format.find('d')])

    grnd_ht = z[0]
    for line in open(specification, 'r'):
        if "Ground Height" in line:
            grnd_ht = float(line[18:])
            break

    if alt_max is None:
        alt_max = z[-1]
    else:
        alt_max = float(alt_max)

    ht_mask = np.logical_and(grnd_ht <= z, z <= alt_max)

    f, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 4]}, figsize=(12, 5))
    
    ax[0].grid(color='k', linestyle='--', linewidth=0.5)
    ax[1].grid(color='k', linestyle='--', linewidth=0.5)

    ax[0].set_ylim(grnd_ht, alt_max)
    ax[0].set_ylabel("Altitude [km]")
    ax[0].set_xlabel("Sound Speed [m/s]")

    ax[1].set_ylim(grnd_ht, alt_max)
    ax[1].set_xlabel("Wind Speed [m/s]")

    ax[1].yaxis.set_ticklabels([])

    ax[2].yaxis.set_label_position("right")
    ax[2].yaxis.tick_right()

    ax[2].set_xlim(-180.0, 180.0)
    ax[2].set_ylim(0.0, 50.0)
    ax[2].set_xticks((-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0, 180.0))
    ax[2].set_xticklabels(["S", "SW", "W", "NW", "N", "NE", "E", "SE", "S"])
    ax[2].set_xlabel("Propagation Direction")
    ax[2].set_ylabel("Inclination [deg]")

    ax[0].plot(c[ht_mask], z[ht_mask], '-k', linewidth=3.0)
    ax[1].plot(u[ht_mask], z[ht_mask], '-b', linewidth=3.0, label='Zonal')
    ax[1].plot(v[ht_mask], z[ht_mask], '-r', linewidth=3.0, label='Merid.')
    ax[1].legend(fontsize='small')

    incl_vals = np.arange(0.0, 50.0, 0.2)
    for az in np.arange(-180.0, 180.0, 1.0):
        ceff = (c + u * np.sin(np.radians(az)) + v * np.cos(np.radians(az)))[ht_mask]
        refract_ht = [z[ht_mask][np.min(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0])] if len(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0]) > 0 else alt_max for incl in incl_vals]
        sc = ax[2].scatter([az] * len(refract_ht), incl_vals, c=refract_ht, cmap=cm.jet_r, marker="s", s=5.0, alpha=0.75, edgecolor='none', vmin=grnd_ht, vmax=120.0)

    f.colorbar(sc, ax=[ax[2]], location='top', label="Estimated Refraction Altitude [km]")

    plt.show()

