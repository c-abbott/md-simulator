"""
CMod Project: Velocity Verlet time integration of
a 3-D system of Argon particles subject to a Leonard-Jones 
potential.

Requires inputs consisting of:
1.) A parameter file
2.) A file to write the data for the MSD plot
3.) A file to write the data for the RDF plot

Produces outputs consisting of:
1.) A plot of the kinetic, potential and total energy of the system 
2.) A plot of the mean squared displacement of the particles
3.) A plot of the radial distribution function
4.) 4 seperate files containing the raw data for the RDF, MSD, VMD and energy plot

Authors: T.E.Bruggi, C.B.Abbott
Version: 11 Mar 2019
"""
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from Particle3D import Particle3D
import pbc
import mdutilities as MDU
import energyforce as ef
import obs

# four functions used to write all data to files


def vmd_writer(output_file, count, particles):
    output_file.write("{0:d}\nPoint = {1}\n".format(len(particles), count))
    # provide applicable format for VMD
    for j in range(len(particles)):
        output_file.write('{0!s}\n'.format(particles[j].__str__()))


def energy_writer(output_file, time, KE_list, PE_list, tot_energy_list):
    output_file.write("#Time, KE, PE, Total Energy\n")

    for i in range(len(tot_energy_list)):
        # round to 3 sig. figures to correctly format time outputs
        output_file.write("%.3f" % time[i] + ", " + "%.7f" % KE_list[i] + ", " + "%.7f" % PE_list[i] + ", " + "%.7f" % tot_energy_list[i] + "\n")
        output_file.close


def msd_writer(output_file, time, MSD_list):
    output_file.write("#Time, MSD\n")

    for i in range(len(MSD_list)):
        # provide a consistent format for outputs
        output_file.write("%.3f" % time[i] + ", " + "%.7f" % MSD_list[i] + "\n")
        output_file.close


def rdf_writer(output_file, bins_data, counts_data):
    output_file.write("#bin, rdf\n")

    for i in range(len(counts_data)):
        output_file.write("%.3f" % bins_data[i] + ", " + str(counts_data[i]) + "\n")
        output_file.close


def main():
    # demand inputs in order for our code to produce outputs
    if len(sys.argv) != 6:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <parameters file>" + "<vmd file>" + "<rdf file>" + "<msd file>" + "<energy file>")
        quit()
    else:
        infile_parameters = sys.argv[1]
        outfile_VMD = sys.argv[2]
        outfile_RDF = sys.argv[3]
        outfile_MSD = sys.argv[4]
        outfile_energy = sys.argv[5]

    # define filehandles
    output_file_VMD = open(outfile_VMD, "w")
    output_file_RDF = open(outfile_RDF, "w")
    output_file_MSD = open(outfile_MSD, "w")
    output_file_energy = open(outfile_energy, "w")

    # open input file and assinging parameters
    with open(infile_parameters, "r") as input_file:
        # read the lines of the input data file
        line = input_file.readline()
        items = line.split(", ")

        natoms = int(items[0])  # no. of particles
        steps = int(items[1])  # no. of simulation steps
        dt = float(items[2])  # timestep
        temp = float(items[3])  # reduced temperature
        rho = float(items[4])  # reduced density
        rc = float(items[5])  # LJ cutoff distance
        n = int(items[6])  # nth timestep for data collection
        nbins = int(items[7])  # number of hist bins for rdf plot

    # create a list to hold all particle objects
    particles = []
    for i in range(natoms):
        particles.append(Particle3D("Argon" + str(i), 1.0, np.array([0, 0, 0]), np.array([1, 1, 1])))

    # set initial positions and velocities of particles in our simulation
    cell = MDU.set_initial_positions(rho, particles)
    MDU.set_initial_velocities(temp, particles)

    # storing initial positions for MSD
    for i in range(natoms):
        particles[i].ini_position = particles[i].position

    # cell_length used for period boundary conditions
    cell_length = cell[0]

    # initialise times for iterative purposes
    time = 0.0
    count = 0

    # initialise data lists for plotting
    time_list = []
    tot_energy_list = []
    KE_list = []
    PE_list = []
    MSD_list = []
    RDF_list = []           # create empty list of separations

    # initialise force values
    forces = ef.force_sum(particles, rc, cell_length)

    # write initial values to VMD file
    vmd_writer(output_file_VMD, count, particles)

    # start the time integration loop
    for i in range(steps):

        # update particle positions
        for j in range(natoms):
            particles[j].leap_position(dt, forces[:, j])
            #pbc.pbc_correction(particles[j].position, cell_length)

        # calculate new forces
        forces_new = ef.force_sum(particles, rc, cell_length)

        # update particle velocities
        for j in range(natoms):
            particles[j].leap_velocity(dt, 0.50 * (forces_new[:, j] + forces[:, j]))
        forces = forces_new

        # append observable data every nth timestep (MSD, RDF and energy)
        if count % n == 0:
            MSD_list.append(OBS.get_MSD(particles, cell_length))
            OBS.get_RDF(particles, cell_length, RDF_list)
            total_energy, KE, PE = ef.get_energy(particles, rc, cell_length)
            tot_energy_list.append(total_energy)
            KE_list.append(KE)
            PE_list.append(PE)
            time_list.append(time)

        # write VMD data to file every timestep
        vmd_writer(output_file_VMD, count, particles)

        # increase time
        time += dt
        count += 1

    # Post-simulation:

    # caluclate the relative energy fluctuation of the system
    e_uncertainty = abs((max(tot_energy_list) - min(tot_energy_list)) / (sum(tot_energy_list) / len(tot_energy_list)))
    e_uncertainty = round(e_uncertainty, 7)
    print("The relative total energy fluctiation of the simulation is: " + str(e_uncertainty))

    # determine the diffusion constant
    D = (MSD_list[-1] - MSD_list[0]) / (6 * time)
    print("The diffusion constant (in reduced units) is: " + str(round(D, 7)))

    # estimate the equilibration time
    KE_average = sum(KE_list) / len(KE_list)
    t0 = 0
    for i in range(len(KE_list)):
        if KE_list[i] > KE_average:
            t0 += dt * n
        else:
            break
    print("The equilibration time (in reduced units) of the system is: " + str(round(t0, 5)))

    # plot kinetic, potential and total energy
    plt.title('Total, Kinetic and Potential Energy of the system vs Time')
    plt.xlabel('Time (reduced units)')
    plt.ylabel('Energy (reduced units)')
    plt.plot(time_list, tot_energy_list, "k-", label="Total Energy", linewidth=2)
    plt.plot(time_list, KE_list, "b-", label="Kinetic Energy")
    plt.plot(time_list, PE_list, "r-", label="Potential Energy")
    plt.legend()
    plt.show()

    # plot MSD
    plt.title('Mean Squared Displacement')
    plt.xlabel('Time (reduced units)')
    plt.ylabel('MSD')
    plt.plot(time_list, MSD_list)
    plt.show()

    # extract RDF data and the bins
    RDF_data, bins_data, bars = plt.hist(RDF_list, nbins)
    plt.close()

    # normalize our RDF (call a function rather than letting plt.hist do it)
    RDF_norm = OBS.normalisation(bins_data, RDF_data, RDF_list)

    # plot RDF
    plt.title('Radial Distribution Function')
    plt.xlabel('Separation (reduced units)')
    plt.ylabel('g(r) - (RDF)')
    plt.plot(bins_data[:-1], RDF_norm, color="k")
    plt.xlim(0), plt.ylim(0)
    plt.show()

    # write to output files
    msd_writer(output_file_MSD, time_list, MSD_list)
    energy_writer(output_file_energy, time_list, KE_list, PE_list, tot_energy_list)
    rdf_writer(output_file_RDF, bins_data[:-1], RDF_norm)


main()
