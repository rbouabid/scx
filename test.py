import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import root_numpy as rnp

import common

file_name = '/home/ryan/Documents/PionChargeExchangeMerged.root'
tree_name = 'PionChargeExchange/tpc'
branch_list = [
    'energy_deposited',
    'number_photon_conversions',
    'photon_start_energy',
    'photon_convert',
    'shower_tree_energy',
    'shower_tree_containment',
    'photon_start_x',
    'photon_start_y',
    'photon_start_z',
    'photon_end_x',  
    'photon_end_y',
    'photon_end_z',
    'primary_start_pe',
    'neutral_pion_start_pe',
]

tree_array = rnp.root2array(
    file_name,
    tree_name,
    branch_list,
    )

number_entries = tree_array.shape[0]

# Accessing Tree Variables and Decs
energy_deposited_array = tree_array['energy_deposited']
number_photon_conversions = tree_array['number_photon_conversions']
gamma_energy_array = tree_array['photon_start_energy']
shower_tree_energy_array = tree_array['shower_tree_energy']
shower_tree_containment_array = tree_array['shower_tree_containment']
photon_start_x_array = tree_array['photon_start_x']
photon_dist_to_edge_array = []
photon_angle_array = []
photon_cos_array = []
pion_mass_array = []
mc_pion_mass_array = []
primary_energy_array = []
neutral_pion_energy_array = []
primary_momentum_array = []
mc_neutral_pion_momentum_array = []
neutral_pion_momentum_array = []

flag = number_photon_conversions > 1

for entry in tree_array:
    photon_0_start = np.array([
        entry['photon_start_x'][0],
        entry['photon_start_y'][0],
        entry['photon_start_z'][0]])
    photon_0_end = np.array([
        entry['photon_end_x'][0],
        entry['photon_end_y'][0],
        entry['photon_end_z'][0]])
    photon_1_start = np.array([
        entry['photon_start_x'][1],
        entry['photon_start_y'][1],
        entry['photon_start_z'][1]])
    photon_1_end = np.array([
        entry['photon_end_x'][1],
        entry['photon_end_y'][1],
        entry['photon_end_z'][1]])

    photon_0_vector = photon_0_end - photon_0_start
    photon_1_vector = photon_1_end - photon_1_start
    photon_angle = common.angle_between(photon_0_vector, photon_1_vector)
    photon_cosine = common.cosine_between(photon_0_vector, photon_1_vector)
    
    shower_0_energy = entry['shower_tree_energy'][0]
    shower_1_energy = entry['shower_tree_energy'][1]

    photon_0_mc_energy = entry['photon_start_energy'][0] 
    photon_1_mc_energy = entry['photon_start_energy'][1] 
   
    primary_energy = entry['primary_start_pe'][3]
    primary_momentum = np.sqrt((entry['primary_start_pe'][0])*(entry['primary_start_pe'][0]) + (entry['primary_start_pe'][1])*(entry['primary_start_pe'][1]) + (entry['primary_start_pe'][2])*(entry['primary_start_pe'][2]))
    neutral_pion_energy = entry['neutral_pion_start_pe'][3]
    mc_neutral_pion_momentum = np.sqrt((entry['neutral_pion_start_pe'][0])*(entry['neutral_pion_start_pe'][0]) + (entry['neutral_pion_start_pe'][1])*(entry['neutral_pion_start_pe'][1]) + (entry['neutral_pion_start_pe'][2])*(entry['neutral_pion_start_pe'][2]))
       
 
    photon_dist_to_edge_array.append([
        common.distance_from_tpc_edge(photon_0_start, photon_0_end),
        common.distance_from_tpc_edge(photon_1_start, photon_1_end),
        ])
    photon_cos_array.append(photon_cosine)
    photon_angle_array.append(photon_angle)
    pion_mass_array.append(
        common.pion_mass(shower_0_energy, shower_1_energy, photon_cosine))
    mc_pion_mass_array.append(
        common.pion_mass(photon_0_mc_energy, photon_1_mc_energy, photon_cosine))
    neutral_pion_momentum_array.append(
        common.neutral_pion_momentum(shower_0_energy, shower_1_energy, photon_angle))

    primary_energy_array.append(primary_energy)
    primary_momentum_array.append(primary_momentum)
    neutral_pion_energy_array.append(neutral_pion_energy)
    mc_neutral_pion_momentum_array.append(mc_neutral_pion_momentum)

photon_dist_to_edge_array = np.array(photon_dist_to_edge_array)
photon_angle_array = np.array(photon_angle_array)
photon_cos_array = np.array(photon_cos_array)
pion_mass_array = np.array(pion_mass_array)
mc_pion_mass_array = np.array(mc_pion_mass_array)
primary_energy_array = np.array(primary_energy_array)
primary_momentum_array = np.array(primary_momentum_array)
neutral_pion_energy_array = np.array(neutral_pion_energy_array)
neutral_pion_momentum_array = np.array(neutral_pion_momentum_array)
mc_neutral_pion_momentum_array = np.array(mc_neutral_pion_momentum_array)

#print photon_dist_to_edge_array[flag]

#print number_photon_conversions
#print flag
#
#print gamma_energy_array[flag]
#print shower_tree_energy_array[flag]
#print shower_tree_containment_array[flag]
#
#print np.min(shower_tree_containment_array[flag], axis=1)

#Plotting

#Gamma Count
f, ax = plt.subplots()
ax.hist(number_photon_conversions, bins=3, range=(0, 2), ec = 'none', color = 'g')
ax.set_xlabel('Number of photons converted in TPC')
ax.set_ylabel('Number of Entries')
ax.set_title('Gamma Count')
ax.legend()
plt.show()

#Primary Particle (+- pions) start energy
f, ax = plt.subplots()
ax.hist(primary_energy_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
ax.set_xlabel('Primary particle start PE(MeV)')
ax.set_ylabel('Number of Entries')
ax.set_title('Primary startE')
plt.show()

#Primary particle start momentum (magnitude)
f, ax = plt.subplots()
ax.hist(primary_momentum_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
ax.set_xlabel('Primary particle start Momentum(magnitude in MeV)')
ax.set_ylabel('Number of Entries')
ax.set_title('Primary startP')
plt.show()

#Pi0 start energy
f, ax = plt.subplots()
ax.hist(neutral_pion_energy_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
ax.set_xlabel('Neutral pion start energy(MeV)')
ax.set_ylabel('Number of Entries')
ax.set_title('Pi0 startE')
plt.show()

#"Truth" Pi0 start momentum (magnitude)
f, ax = plt.subplots()
ax.hist(mc_neutral_pion_momentum_array[flag], bins=400, range=(0,4000), ec = 'none', color = 'g', alpha = 0.5, label = 'truth')
ax.hist(neutral_pion_momentum_array[flag], bins=400, range=(0,4000), ec = 'none', color = 'r', alpha = 0.5, label = 'deposited')
ax.set_xlabel('Neutral pion start Momentum(magnitude)(MeV)')
ax.set_ylabel('Number of Entries')
ax.set_title('Pi0 momentum')
ax.legend(bbox_to_anchor=(0.9, 0.9))
plt.show()

#MC shower energy

f, ax = plt.subplots()
ax.hist(shower_tree_energy_array[flag].flatten(), bins=75, range=(0, 750), ec = 'none', color = 'g')
ax.set_xlabel('gamma energy deposited(MeV)')
ax.set_ylabel('Number of Entries')
ax.set_title('Energy Deposition')
plt.show()

#Containment
f, ax = plt.subplots()
ax.hist(shower_tree_containment_array[flag].flatten(), bins=10, range=(0, 1), ec = 'none', color = 'g')
ax.set_xlabel('Containment')
ax.set_ylabel('Number of Entries')
ax.set_title('Containment Ratio')
plt.show()

#Forward dist
f, ax = plt.subplots()
ax.hist(photon_dist_to_edge_array[flag].flatten(), bins=50, range=(0, 110), ec = 'none', color = 'g')
ax.set_xlabel('Distance(cm)')
ax.set_ylabel('Number of Entries')
ax.set_title('Distance to TPC edge')
plt.show()

#Mass
f, ax = plt.subplots()
plt.hist(pion_mass_array[flag], bins=70, range=(0, 160), ec = 'none', color = 'g', alpha = 0.5, label = 'deposited')
plt.hist(mc_pion_mass_array[flag], bins=70, range=(0, 160), ec = 'none', color = 'r', alpha = 0.5, label = 'truth')
plt.xlabel('Mass(MeV)')
plt.ylabel('Number of Entries')
plt.title('Pion Mass')
ax.legend(bbox_to_anchor=(0.9, 0.9))
plt.show()

#dist vs contain
f, ax = plt.subplots()
h2d, xbins, ybins, im = ax.hist2d(photon_dist_to_edge_array[flag].flatten(), shower_tree_containment_array[flag].flatten(), bins = (55, 20), range = [(0,110),(0,1)], norm=LogNorm())
ax.set_xlabel('Distance to edge(cm)')
ax.set_ylabel('Containment ratio')
ax.set_title('Distance vs Containment')
f.colorbar(im)
plt.show()

#dist vs mass
f, ax = plt.subplots()
h2d, xbins, ybins, im = ax.hist2d(np.amin(photon_dist_to_edge_array[flag], axis=1), pion_mass_array[flag], bins = (55, 70), range = [(0,110),(0,140)], norm=LogNorm())
ax.set_xlabel('Min distance to edge(cm)')
ax.set_ylabel('Pi0 mass(MeV)')
ax.set_title('minDistance vs Mass')
f.colorbar(im)
plt.show()

#contain vs mass
f, ax = plt.subplots()
h2d, xbins, ybins, im = ax.hist2d(np.amin(shower_tree_containment_array[flag], axis=1), pion_mass_array[flag], bins = (20, 70), range = [(0,1),(0,140)])
ax.set_xlabel('Min containment ratio')
ax.set_ylabel('Pi0 mass(MeV)')
ax.set_title('minContainment vs Mass')
f.colorbar(im)
plt.show()
