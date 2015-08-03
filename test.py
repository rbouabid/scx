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
    'neutral_pion_start_pe'
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
photon_cos_array = []
pion_mass_array = []
primary_energy_array = []
neutral_pion_energy_array = []
primary_momentum_array = []
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
    photon_cosine = common.cosine_between(photon_0_vector, photon_1_vector)
    shower_0_energy = entry['shower_tree_energy'][0]
    shower_1_energy = entry['shower_tree_energy'][1]
    primary_energy = entry['primary_start_pe'][3]
    primary_momentum = np.sqrt((entry['primary_start_pe'][0])*(entry['primary_start_pe'][0]) + (entry['primary_start_pe'][1])*(entry['primary_start_pe'][1]) + (entry['primary_start_pe'][2])*(entry['primary_start_pe'][2]))
    neutral_pion_energy = entry['neutral_pion_start_pe'][3]
    neutral_pion_momentum = np.sqrt((entry['neutral_pion_start_pe'][0])*(entry['neutral_pion_start_pe'][0]) + (entry['neutral_pion_start_pe'][1])*(entry['neutral_pion_start_pe'][1]) + (entry['neutral_pion_start_pe'][2])*(entry['neutral_pion_start_pe'][2]))
       
 
    photon_dist_to_edge_array.append([
        common.distance_from_tpc_edge(photon_0_start, photon_0_end),
        common.distance_from_tpc_edge(photon_1_start, photon_1_end),
        ])
    photon_cos_array.append(photon_cosine)
    pion_mass_array.append(
        common.pion_mass(shower_0_energy, shower_1_energy, photon_cosine))
    primary_energy_array.append(primary_energy)
    primary_momentum_array.append(primary_momentum)
    neutral_pion_energy_array.append(neutral_pion_energy)
    neutral_pion_momentum_array.append(neutral_pion_momentum)

photon_dist_to_edge_array = np.array(photon_dist_to_edge_array)
photon_cos_array = np.array(photon_cos_array)
pion_mass_array = np.array(pion_mass_array)
primary_energy_array = np.array(primary_energy_array)
primary_momentum_array = np.array(primary_momentum_array)
neutral_pion_energy_array = np.array(neutral_pion_energy_array)
neutral_pion_momentum_array = np.array(neutral_pion_momentum_array)

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
n, bins, patches = plt.hist(number_photon_conversions, bins=3, range=(0, 2), ec = 'none', color = 'g')
plt.xlabel('Number of photons converted in TPC')
plt.ylabel('Number of Entries')
plt.title('Gamma Count')
plt.show()

#Primary Particle (+- pions) start energy
n, bins, patches = plt.hist(primary_energy_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
plt.xlabel('Primary particle start PE(MeV)')
plt.ylabel('Number of Entries')
plt.title('Primary startE')
plt.show()

#Primary particle start momentum (magnitude)
n, bins, patches = plt.hist(primary_momentum_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
plt.xlabel('Primary particle start Momentum(magnitude in MeV)')
plt.ylabel('Number of Entries')
plt.title('Primary startP')
plt.show()

#Pi0 start energy
n, bins, patches = plt.hist(neutral_pion_energy_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
plt.xlabel('Neutral pion start energy(MeV)')
plt.ylabel('Number of Entries')
plt.title('Pi0 startE')
plt.show()

#Pi0 start momentum (magnitude)
n, bins, patches = plt.hist(neutral_pion_momentum_array, bins=400, range=(0,4000), ec = 'none', color = 'g')
plt.xlabel('Neutral pion start Momentum(magnitude)(MeV)')
plt.ylabel('Number of Entries')
plt.title('Pi0 startP')
plt.show()
sys.exit()

n, bins, patches = plt.hist(shower_tree_energy_array[flag].flatten(), bins=75, range=(0, 750), ec = 'none', color = 'g')
plt.xlabel('gamma energy deposited(MeV)')
plt.ylabel('Number of Entries')
plt.title('Energy Deposition')
plt.show()


n, bins, patches = plt.hist(shower_tree_containment_array[flag].flatten(), bins=10, range=(0, 1), ec = 'none', color = 'g')
plt.xlabel('Containment')
plt.ylabel('Number of Entries')
plt.title('Containment Ratio')
plt.show()

n, bins, patches = plt.hist(photon_dist_to_edge_array[flag].flatten(), bins=50, range=(0, 110), ec = 'none', color = 'g')
plt.xlabel('Distance(cm)')
plt.ylabel('Number of Entries')
plt.title('Distance to TPC edge')
plt.show()

n, bins, patches = plt.hist(pion_mass_array[flag], bins=70, range=(0, 140), ec = 'none', color = 'g')
plt.xlabel('Mass(MeV)')
plt.ylabel('Number of Entries')
plt.title('Pion Mass')
plt.show()

plt.hist2d(photon_dist_to_edge_array[flag].flatten(), shower_tree_containment_array[flag].flatten(), bins = (55, 20), range = [(0,110),(0,1)], norm=LogNorm())
plt.xlabel('Distance to edge(cm)')
plt.ylabel('Containment ratio')
plt.title('Distance vs Containment')
plt.colorbar()
plt.show()

plt.hist2d(np.amin(photon_dist_to_edge_array[flag], axis=1), pion_mass_array[flag], bins = (55, 70), range = [(0,110),(0,140)], norm=LogNorm())
plt.xlabel('Min distance to edge(cm)')
plt.ylabel('Pi0 mass(MeV)')
plt.title('minDistance vs Mass')
plt.colorbar()
plt.show()

plt.hist2d(np.amin(shower_tree_containment_array[flag], axis=1), pion_mass_array[flag], bins = (20, 70), range = [(0,1),(0,140)])
plt.xlabel('Min containment ratio')
plt.ylabel('Pi0 mass(MeV)')
plt.title('minContainment vs Mass')
plt.colorbar()
plt.show()
