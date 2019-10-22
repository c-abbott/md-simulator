"""
A module used to fix particle seperations using 
MIC and particle positions using PBC in an MD simulation.

Authors: C.B.Abbott T.E.Bruggi
Version: 03/2019
"""

import math
import numpy as np


def pbc_correction(position_vector, cell_length):

    # evaluating the remainder for each component of the postion vector
    # using the modulo operator to give me the corrected position of the particle
    # within the original cell (PBC)
    position_vector = np.mod(position_vector, cell_length)

    return position_vector


def mic_correction(separation, cell_length):

    # ensuring minimum image convention is obeyed for all particle separations

    separation = np.mod(separation + cell_length / 2, cell_length) - cell_length / 2
    return separation
