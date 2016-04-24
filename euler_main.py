#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy as np
from initialize import initial_setup
import fluxes as flux
from boundary_conditions import apply_boundary_conditions
from timestepping import set_timestep
from smoothing import smooth
from convergence import check_convergence
def main():

    # Setup the grid and compute the initial flow solution
    primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters = initial_setup()


    # Set the timestepping
    set_timestep()

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
    while continue_flag == 0:




        set_other_flow_properties()

        apply_boundary_conditions()


        # Set the fluxes!
        flux.set_fluxes(primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters)

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
