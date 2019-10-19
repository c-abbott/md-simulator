"""
A class for a particle paramterized in 3D

Authors: C.B.Abbott T.E.Bruggi
Version: Mar 2019
"""
import numpy as np
import pbc


class Particle3D(object):
    """
    Class to describe 3D particles

    Properties:
    mass(float) - particle mass
    label(string) - distinguishing label
    postion(NumPy array) - postiion in 3D space
    veloctiy(NumPy array) - velocity in 3D space

    Methods:
    * postion output for each particle
    * kinetic energy
    * first order velocity update
    * first and second order postion updates
    * particle creation from file entry
    * vector seperation of two particles
    """

    def __init__(self, label, mass, pos, vel):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param mass: mass as float
        :param pos: position as array
        :param vel: velocity as array
        """
        self.label = str(label)
        self.mass = float(mass)
        if self.mass < 0:
            print("Error: Check mass input. Mass cannot be negative.")
            quit()
        self.ini_position = pos 
        self.position = pos
        self.velocity = vel
       
    
    def __str__(self):
        """
        String method that prints the particle in the following format:
        <label> <pos_x> <pos_y> <pos_z> 
        """
        
        return (str(self.label) + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2]))
        
    def kinetic_energy(self):
        """
        Return a particle's kinetic energy as
        1/2*mass*velocity^2
        """
        
        return 0.5*self.mass*(np.inner(self.velocity, self.velocity))
    
    # time integration methods

    def leap_velocity(self, dt, force):
        """
        First order velocity update,
        v(t+dt) = v(t) + (F(t)*dt)/m

        :param dt: timestep as float
        :param force: force on particle as NumPy array
        """
        self.velocity += dt*force/self.mass


    def leap_position(self, dt, force):
        """
        Second order position update
        r(t+dt) = r(t) + dt*v(t) + dt^2*F(t)/2m

        :param dt: timestep as float
        :param force: force on particle as a NumPy array
        """
        # do not use += to avoid having inipositions redefined in MSD
        self.position = self.position + dt*self.velocity + ((dt)**2)*force/(2*self.mass)

    @staticmethod
    def particle_separation(p1, p2, cell_length):
        """
        Class method:
        Calculates the vector separation of a particle (p1) with respect to 
        a particle (p2) in 3D as a NumPy array, corrected to minimum 
        image convention.
        
        :param p1: object particle 1
        :param p2: object particle 2
        
        :return: position of p1 - position of p2
        """
        separation = p1.position - p2.position
        rel_separation = pbc.mic_correction(separation, cell_length)
        
        return(rel_separation)
        
