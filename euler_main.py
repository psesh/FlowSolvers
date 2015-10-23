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

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set_fluxes()

Quantities to compute:
~~~~~~~~~~~~~~~~~~~~~~~
fluxi_mass(i,j) |  fluxj_mass(i,j)
fluxi_xmom(i,j) |  fluxj_xmom(i,j)
fluxi_ymom(i,j) |  fluxj_ymom(i,j)
fluxi_enth(i,j) |  fluxj_enth(i,j)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
def set_fluxes():

    #--------------------------------------------------------------------------
    # Calculate the mass fluxes
    #--------------------------------------------------------------------------
    # Mass flux across each "i" face of elements.
    # Also compute total mass flow rate across each "i" line
    for i in range(0, ni):
        flow(i) = 0.0
        for j in range(0, nj - 1):
            global fluxi_mass(i,j) = 0.5 *( (rovx(i,j) + rovx(i, j+1) ) * dlix(i,j) +
                             (rovy(i,j) + rovy(i, j+1) )* dliy(i,j) )
            flow(i) = flow(i) + fluxi_mass(i,j)

    # Now the mass flux across each "j" face
    for i in range(0, ni - 1):
        for j in range(1, nj):
            global fluxj_mass(i,j) = 0.5 *( (rovx(i,j) + rovx(i + 1, j)) * dljx(i,j) +
                            (rovy(i,j) + rovy(i + 1, j)) * dljy(i,j) )


    # Set the mass fluxes through each j={1,...,nj} to be 0 as these are solid
    # surfaces. Not necessary to resolve the velocity parallel to these surfaces
    for i in range(0, ni - 1):
        global fluxj_mass(i,1) = 0.0
        global fluxj_mass(i,nj) = 0.0

    #--------------------------------------------------------------------------
    # Calculate the fluxes of X-momentum
    #--------------------------------------------------------------------------
    # in the "i" direction
    for i in range(0, ni):
        for j in range(0,nj-1):
            global fluxi_xmom(i,j) = 0.5 * ( fluxi_mass(i,j)*(vx(i,j) + vx(i, j+1) ) +
                              (p(i,j) + p(i, j + 1)) * dlix(i,j) )
    # in the "j" direction
    for i in range(0, ni - 1):
        for j in range(0, nj):
            global fluxj_xmom(i,j) = 0.5 * ( fluxj_mass(i,j)*(vx(i,j) + vx(i + 1, j)) +
                             (p(i,j) + p(i+1,j))*dljx(i,j) )

    #--------------------------------------------------------------------------
    # Calculate the fluxes of Y-momentum
    #--------------------------------------------------------------------------
    # in the "i" direction
    for i in range(0, ni):
        for j in range(0, nj - 1):
            fluxi_ymom(i,j) = 0.5 * ( fluxi_mass(i,j) * (vy(i,j) + vy(i, j+1)) +
                                (p(i,j) + p(i,j+1))*dliy(i,j) )

    # in the "j" direction
    for i in range(0, ni - 1):
        for j in range(0, nj):
            fluxj_ymom(i,j) = 0.5 * (fluxj_mass(i,j) * (vy(i,j) + vy(i+1,j)) +
                                (p(i,j) + p(i+1, j)) * dljy(i,j))

    #--------------------------------------------------------------------------
    # Calculate the fluxes of enthalpy
    #--------------------------------------------------------------------------
    # in the "i" direction
    for i in range(0, ni):
        for j in range(0, nj - 1):
            fluxi_enth(i,j) = 0.5 * (fluxi_mass(i,j) + (hstag(i,j) + hstag(i,j+1)))

    # in the "j" direction
    for i in range(0, ni - 1):
        for j in range(0, nj):
            fluxj_enth(i,j) = 0.5 * (fluxj_mass(i,j) + (hstage(i,j) + hstag(i+1, j)))

"""
subroutine sums the fluxes for each element, calculates the changes in these
variable "prop" and distributes them to the four corners of the element
"""
def sum_fluxes(iflux, jflux, prop, delprop):
