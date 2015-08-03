#!/usr/bin/env python

from __future__ import division
import sys
import argparse

import numpy as np
from scipy.optimize import curve_fit  # for fitting
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.colors import LogNorm
import root_numpy as rnp

import common

parser = argparse.ArgumentParser(description="Plot from ROOT file.")
parser.add_argument("file", type=str, help="path to ROOT file")
args = parser.parse_args()

file_path = args.file

print("Plotting from {}".format(file_path))

branch_list = [
    'neutral_pion_start_pe',
    'number_photon_conversions',
    'photon_start_x',
    'photon_start_y',
    'photon_start_z',
    'photon_end_x',
    'photon_end_y',
    'photon_end_z',
    'photon_start_px',
    'photon_start_py',
    'photon_start_pz',
    'photon_start_energy',
    'photon_convert',
    'shower_tree_energy',
    'shower_tree_containment',
    ]

arr = rnp.root2array(file_path, 'PionChargeExchange/tpc', branch_list)

number_photon_conversions = arr['number_photon_conversions']

flag = number_photon_conversions == 2

neutral_pion_start_pe = arr['neutral_pion_start_pe'][flag]
photon_start_x = arr['photon_start_x'][flag]
photon_start_y = arr['photon_start_y'][flag]
photon_start_z = arr['photon_start_z'][flag]
photon_end_x = arr['photon_end_x'][flag]
photon_end_y = arr['photon_end_y'][flag]
photon_end_z = arr['photon_end_z'][flag]
photon_start_px = arr['photon_start_px'][flag]
photon_start_py = arr['photon_start_py'][flag]
photon_start_pz = arr['photon_start_pz'][flag]
photon_start_energy = arr['photon_start_energy'][flag]
photon_convert = arr['photon_convert'][flag]
shower_tree_energy = arr['shower_tree_energy'][flag]

photon_vector_x = photon_end_x - photon_start_x
photon_vector_y = photon_end_y - photon_start_y
photon_vector_z = photon_end_z - photon_start_z

neutral_pion_mc_momentum = []
neutral_pion_momentum = []
photon_mc_angle = []
photon_mc_cosine_angle = []

number_events = number_photon_conversions[flag].size

for i in xrange(number_events):

    # get neutral pion Monte Carlo truth momentum
    neutral_pion_mc_momentum.append(neutral_pion_start_pe[i][-1])

    # get directional vectors of photons from Monte Carlo truth information
    photon_1_vector = np.array([
        photon_vector_x[i][0], photon_vector_y[i][0], photon_vector_z[i][0] ])
    photon_2_vector = np.array([
        photon_vector_x[i][1], photon_vector_y[i][1], photon_vector_z[i][1] ])

    # get Monte Carlo truth deposited energy of shower
    shower_1_mc_energy = shower_tree_energy[i][0]
    shower_2_mc_energy = shower_tree_energy[i][1]

    # get opening angle between the two photons
    opening_angle = common.angle_between(photon_1_vector, photon_2_vector)
    cosine_angle = common.cosine_between(photon_1_vector, photon_2_vector)

    photon_mc_angle.append(opening_angle * 180.0 / np.pi)
    photon_mc_cosine_angle.append(cosine_angle)

    # get neutral pion momentum
    neutral_pion_momentum.append(
        common.neutral_pion_momentum(
            shower_1_mc_energy, shower_2_mc_energy, opening_angle)
        )

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#h2d, xbins, ybins, im = ax.hist2d(photon_mc_angle,
#                                  neutral_pion_mc_momentum,
#                                  bins=[ 180, 180 ],
#                                  range=[(0.0, 180.0), (0.0, 1800.0)],
#                                  norm=LogNorm())
ax.scatter(photon_mc_angle, neutral_pion_mc_momentum, color='k', s=1)
ax.set_xlim([0, 180])
ax.set_ylim([0, 3000])
ax.set_xlabel(r"${\theta}_{{\gamma}{\gamma}}$")
ax.set_ylabel(r"${\pi}^{0}$ momentum [MeV]")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(photon_mc_cosine_angle, neutral_pion_mc_momentum, color='k', s=1)
ax.set_xlim([-1, 1])
ax.set_ylim([0, 3000])
ax.set_xlabel(r"$\cos\,{\theta}_{{\gamma}{\gamma}}$")
ax.set_ylabel(r"${\pi}^{0}$ momentum [MeV]")
plt.show()
