"""
A module used for caluculating the observables of a system, in particular
the normalised radial distribution function and the mean square
displacement.

The functions consist of:
- MSD calculator
- RDF calculator
- RDF normalisation

Authors: C.B.Abbott and T.E.Bruggi
Version: 03/2019
""" 
import numpy as np
from Particle3D import Particle3D
import pbc

def get_MSD(particles, cell_length):
    msd = 0.0
    for i in range(len(particles)):
        separation = np.linalg.norm(particles[i].position - particles[i].ini_position)

        # apply periodic boundary conditions corrections
        separation = pbc.mic_correction(separation, cell_length)
        msd += (1/len(particles))*separation**2

    return msd

def get_RDF(particles, cell_length, RDF_list):
    count = 0
    # calculating all pairwise seperations whilst
    # being conscious of not overcounting
    for j in range(len(particles)):
        for k in range(j+1, len(particles)):
            # find magnitude of each separation
            rdf = np.linalg.norm(Particle3D.particle_separation(particles[j], particles[k], cell_length))
            RDF_list.append(rdf)
            count += 1

    # return a list of separations
    return RDF_list

def normalisation(bins_data, RDF_data, RDF_list):
    # determine the bin size
    dr = bins_data[1] - bins_data[0]

    count = sum(RDF_data)
    # normalize RDF
    RDF_norm = RDF_data/(count*dr)

    return RDF_norm
