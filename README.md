# md-simulator
Lennard-Jones Simulation For Argon

The code uses a velocity verlet algorithim in order to 
simulate the time evolution of Argon particles confined to
a box using periodic boundary conditions and minimum 
image convention.

It will produce outputs consisting of:
1.) A plot of the kinetic, potential and total energy of the system 
2.) A plot of the mean squared displacement of the particles
3.) A plot of the radial distribution function
4.) A file used to visualise the system in Visual Molecular Dynamics (VMD)
5.) 4 seperate files containing the raw data for the RDF, MSD, VMD and energy plot


In order to run the code you will require 6 entries into
the bash terminal (not including python3) with the usage:
python3 LJArgon.py <parameters file> <vmd file> <rdf file> <msd file> <energy file>

An example of what would be required in the bash terminal would be:
python3 LJArgon.py parameters.txt vmd.xyz rdf_data.txt msd_data.txt energy_data.txt

You can edit the parameters of the simulation by editing the "parameters.txt" file,
I recommend for the three phases of Argon the following temperatures and densities:

SOLID: T = 0.1   rho = 1.0   natoms = 108  LJ cut off distance = 2.5

LIQUID: T = 1.2  rho = 0.67  natoms = 80   LJ cut off distance = 3.0

GAS: T = 1.8  rho = 0.23     natoms = 80   LJ cut off distance = 3.0

All units in the MD simulation are in reduced units.

The other simulation parameters (excluding T and rho) consist of:
no. of particles, no. of simulation steps, timestep,  
LJ cutoff distance, nth timestep, number of histogram bins
which are free to choose as you wish.

Authors: T.E.Bruggi and C.B.Abbott
