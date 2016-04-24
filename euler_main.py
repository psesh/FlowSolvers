#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy as np
from initialize import initial_setup
import fluxes as flux
from boundary_conditions import apply_boundary_conditions
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
    primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters = initial_setup()
    areas = grid_parameters[7] # Will need this for later!

    # Set the time-step
    step = set_timestep(primary_variables, secondary_variables, boundary_conditions, grid_parameters):

    # Time-marching loop
    continue_flag = 0


    display_item =
    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                   TUNGSTEN
                                  ver. 1.0.0

                        Copyright (c) 2016 by Pranay Seshadri
                      University of Cambridge, Cambridge, U.K.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    print(display_item)
    nsteps = 10

    # Two for-loops
    for step_number in range(1, nsteps):
        for nrkut in range(1, 5):

            # setting frkut
            frkut = 1.0/(1.0 + 4 - nrkut)

            # Set "secondary" variables!
            secondary_variables = set_other_variables(primary_variables, secondary_variables, boundary_conditions, grid_parameters):

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

            # Sum the fluxes
            ro, del_ro = sum_fluxes(flux_i_mass, flux_j_mass, ro, ro_start, del_ro, frkut, step, areas)
            ro_vel_x, del_ro_vel_x = sum_fluxes(flux_i_xmom, flux_j_xmom, ro_vel_x, ro_vel_x_start, del_ro_vel_x, frkut, step, areas)
            ro_vel_y, del_ro_vel_y = sum_fluxes(flux_i_ymom, flux_j_ymom, ro_vel_y, ro_vel_y_start, del_ro_vel_y, frkut, step, areas)
            ro_energy, del_ro_energy = sum_fluxes(flux_i_enthalpy, flux_j_enthalpy, ro_energy, ro_energy_start, del_ro_energy, frkut, step, areas)

            # FINISH THE SMOOTHING!


            # Pack everything up!
            primary_variables[0] = ro
            primary_variables[1] = ro_vel_x
            primary_variables[2] = ro_vel_y
            primary_variables[3] = ro_energy

main()
