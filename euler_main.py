#!/usr/bin/python
import numpy as np
from initialize import initial_setup
import fluxes as flux
from boundary_conditions import apply_boundary_conditions, set_other_variables
from timestepping import set_timestep
from smoothing import smooth
from convergence import check_convergence
"""
################################################################################
                                TUNGSTEN FLOW SOLVER

    Simple 2D flow solver to demonstrate the effect of streamline curvature
    in a converging-diverging nozzle by solving the Euler equations on a
    structured mesh. The mesh contours can be varied using a parametric
    definition as given in geometry.py

    Main files:
       1. euler_main.py
       2. boundary_conditions.py
       3. initialize.py
       4. timestepping.py
       5. convergence.py
       6. fluxes.py
       7. geometry.py

################################################################################
"""
def main():


    # Setup the grid and compute the initial flow solution
    nu, nv, nw = 120, 60, 1
    primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters = initial_setup(nu, nv, nw)
    areas = grid_parameters[7] # Will need this for later!

    # Set the time-step
    step = set_timestep(primary_variables, secondary_variables, boundary_conditions, grid_parameters)

    # Output iteration frequency
    niter = 3

    # Initialize some starting values!
    ro_start = np.zeros((nv, nu, nw))
    ro_vel_x_start = np.zeros((nv, nu, nw))
    ro_vel_y_start = np.zeros((nv, nu, nw))
    ro_energy_start = np.zeros((nv, nu, nw))

    # Initialize delta properties
    del_ro = np.zeros((nv, nu))
    del_ro_vel_x = np.zeros((nv, nu))
    del_ro_vel_y = np.zeros((nv, nu))
    del_ro_energy = np.zeros((nv, nu))

    # 'Corrected flow properties'
    corrected_ro = np.zeros((nv, nu, nw))
    corrected_ro_vel_x = np.zeros((nv, nu, nw))
    corrected_ro_vel_y = np.zeros((nv, nu, nw))
    corrected_ro_energy = np.zeros((nv, nu, nw))

    # Display the following message everytime the program is initialized
    display_item = str("""
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                   TUNGSTEN
                                  ver. 1.0.0

                        Copyright (c) 2016 by Pranay Seshadri
                      University of Cambridge, Cambridge, U.K.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~""")
    print(display_item)
    nsteps = 10

    # Two for-loops
    for step_number in range(1, nsteps):

        # Define some starting values
        starting_variables = set_starting_values(primary_variables)

        for nrkut in range(1, 5):

            # Define some starting values!

            # setting frkut
            frkut = 1.0/(1.0 + 4 - nrkut)

            # Set "secondary" variables!
            secondary_variables = set_other_variables(primary_variables, secondary_variables, boundary_conditions, grid_parameters)

            # Enforce boundary conditions
            primary_variables, secondary_variables = apply_boundary_conditions(step_number, primary_variables, secondary_variables, boundary_conditions, grid_parameters)

            # Set the fluxes!
            fluxes = flux.set_fluxes(primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters)

            # Unpack the fluxes
            flux_i_mass = fluxes[0]
            flux_j_mass = fluxes[1]
            flux_i_xmom = fluxes[2]
            flux_j_xmom = fluxes[3]
            flux_i_ymom = fluxes[4]
            flux_j_ymom = fluxes[5]
            flux_i_enthalpy = fluxes[6]
            flux_j_enthalpy = fluxes[7]

            # Unpack the primary variables
            ro = primary_variables[0]
            ro_vel_x = primary_variables[1]
            ro_vel_y = primary_variables[2]
            ro_energy = primary_variables[3]

            # Unpack the starting variables
            ro_start = starting_variables[0]
            ro_vel_x = starting_variables[1]
            

            # Sum the fluxes
            ro, del_ro = flux.sum_fluxes(flux_i_mass, flux_j_mass, ro, ro_start, del_ro, frkut, step, areas)
            ro_vel_x, del_ro_vel_x = flux.sum_fluxes(flux_i_xmom, flux_j_xmom, ro_vel_x, ro_vel_x_start, del_ro_vel_x, frkut, step, areas)
            ro_vel_y, del_ro_vel_y = flux.sum_fluxes(flux_i_ymom, flux_j_ymom, ro_vel_y, ro_vel_y_start, del_ro_vel_y, frkut, step, areas)
            ro_energy, del_ro_energy = flux.sum_fluxes(flux_i_enthalpy, flux_j_enthalpy, ro_energy, ro_energy_start, del_ro_energy, frkut, step, areas)

            # Smoothing!
            ro = smooth(ro, corrected_ro, boundary_conditions, grid_parameters)
            ro_vel_x = smooth(ro_vel_x, corrected_ro_vel_x, boundary_conditions, grid_parameters)
            ro_vel_y = smooth(ro_vel_y, corrected_ro_vel_y, boundary_conditions, grid_parameters)
            ro_energy = smooth(ro_energy, corrected_ro_energy, boundary_conditions, grid_parameters)

            # Pack everything up!
            primary_variables[0] = ro
            primary_variables[1] = ro_vel_x
            primary_variables[2] = ro_vel_y
            primary_variables[3] = ro_energy

        # Write out convergence parameters at every niter iterations
        if(np.mod(step_number, 5) == 0):
            old_primary_variables = check_convergence(primary_variables, old_primary_variables, reference_values, step_number)


main()
