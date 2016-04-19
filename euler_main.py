#!/usr/bin/python
import numpy as np
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

---> need a way to incorporate "common block euler"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def main():

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



" Compute secondary flow variables from the primary ones at each grid point"
" primary variables: ro, rovx, rovy, roe"
" secondary variables: vx, vy, p, hstag"
def set_other_flow_properties():

    for i in range(0, ni):
        for j in range(0, nj):
            global vx(i,j) = rovx(i,j)/ro(i,j)
            global vy(i,j) = rovy(i,j)/ro(i,j)
            e = roe(i,j)/ro(i,j)
            t = (e - 0.5*vx(i,j)**2 - 0.5*vy(i,j)**2)/cv
            global p(i,j) = ro(i,j) * rgas * t
            global hstag(i,j) = cp * tstagin


def set_timestep():

    astag = np.sqrt(gamma * rgas * tstagin)
    umax = astag

    # compute delta_t
    deltat = cfl * dmin/(astag + umax)


"""
This subroutine smooths the variable "prop" (i.e. it adds the artificial viscosity
) by taking (1-SF) * the calculated value of "prop" + SF x (the average of the surrounding
values of "prop"). Here SF is the smoothing factor
"""


