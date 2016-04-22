#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy as np
from initialize import initial_setup
import fluxes as flux
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              "Tungsten"
                              ver. 1.0.0

                    Copyright (c) 2015 by Pranay Seshadri

Features
- 2D finite volume Euler code
- Lax-Friedrichs method
- time derivatives obtained by central differencing
- solution stabilized by artifical viscosity


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def main():
    
    # Setup the grid and compute the initial flow solution
    primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters = initial_setup()
    
    
    # Set the timestepping
    set_timestep()
    
    # Time-marching loop
    continue_flag = 0
    while continue_flag == 0:    
        
        flux.set_fluxes(primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters)
    
    

    # Read in the geometry and flow conditions
    read_data()

    # Generate the grid coordinates, element areas and projected lengths of the
    # sides of the elements
    generate_grid()

    # Check that the areas and the projected lengths are correct
    check_grid()

    # Initial solution
    flow_guess()

    # Set the length of the timestep. Initially this is a constant time step
    # based on a conservative guess of the Mach number
    set_timestep()

    # time stepping for the
    while True:

        set_other_flow_properties()

        apply_boundary_conditions()

        set_fluxes()

        # Sum the fluxes
        sum_fluxes(fluxi_mass, fluxj_mass, ro, delro)
        sum_fluxes(fluxi_enth, fluxj_enth, roe, delroe)
        sum_fluxes(fluxi_xmom, fluxj_xmom, rovx, delrovx)
        sum_fluxes(fluxi_ymom, fluxj_ymom, rovy, delrovy)

        # smoothing
        smooth(ro)
        smooth(rovx)
        smooth(rovy)
        smooth(roe)

        # convergence check here

        if fail_condition:
            break


main()