"""
A module to calculate the total energy and total
force on a particle in an N body simulation.
The potential is modelled as a Lennard-Jones interaction between atoms.

Authors: T.E.Bruggi C.B.Abbott 
Version: 03/2019
"""
import numpy as np
from Particle3D import Particle3D

# function calculating the total energy of the system
# using the potential and kinetic energies
def get_energy(particles, rc, cell_length):
    PE = 0
    KE = 0
    # summing over all particles to calculate
    # all pairwise potentials
    for i in range(len(particles)):
        KE += particles[i].kinetic_energy()
        for j in range(len(particles)):
            if j > i:
                separation = Particle3D.particle_separation(particles[i], particles[j], cell_length)
                r_mag = np.linalg.norm(separation)
                if r_mag > rc:
                    pass
                else:
                    PE += 4.*((1/r_mag)**12 - (1/r_mag)**6) 
            else:
                pass

    total_energy = KE + PE
    return total_energy, KE, PE

# calculating total force due to pairwise
# interactions in the Lennard-Jones Potential
def force_sum(particles, rc, cell_length):
    forces = np.zeros([3, len(particles)])
    for i in range(len(particles)):
        total_force = np.zeros(3)
        for j in range(len(particles)):
            if j != i:
                separation = Particle3D.particle_separation(particles[i], particles[j], cell_length)
                r_mag = np.linalg.norm(separation)
                if r_mag > rc:
                    pass
                else:
                    # calculate the force on particle i due to particle j
                    total_force += 48*((1/r_mag)**14 - 1/2*(1/r_mag)**8)*separation
        forces[:,i] = total_force                  

    return forces

