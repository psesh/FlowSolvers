#!/usr/bin/python
import numpy as np

"""

Class for setting the fluxes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set_fluxes()

Quantities to compute:
~~~~~~~~~~~~~~~~~~~~~~~
fluxi_mass(i,j) |  fluxj_mass(i,j)
fluxi_xmom(i,j) |  fluxj_xmom(i,j)
fluxi_ymom(i,j) |  fluxj_ymom(i,j)
fluxi_enth(i,j) |  fluxj_enth(i,j)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
def set_fluxes(primary_variables, secondary_variables, fluxes, boundary_conditions, grid_parameters):


    # "Unpack" the fluxes
    flux_i_mass = fluxes[0]
    flux_j_mass = fluxes[1]
    flux_i_xmom = fluxes[2]
    flux_j_xmom = fluxes[3]
    flux_i_ymom = fluxes[4]
    flux_j_ymom = fluxes[5]
    flux_i_enthalpy = fluxes[6]
    flux_j_enthalpy = fluxes[7]
    flow = fluxes[8]

    # "Unpack" the primary flow variables
    ro_vel_x = primary_variables[1]
    ro_vel_y = primary_variables[2]

    # "Unpack" the secondary flow variables
    vel_x = secondary_variables[0]
    vel_y = secondary_variables[1]
    pressure = secondary_variables[2]
    enthalpy_stag = secondary_variables[3]

    # "Unpack" the grid parameters
    point_x = grid_parameters[0]
    dlix = grid_parameters[8]
    dliy = grid_parameters[9]
    dljx = grid_parameters[10]
    dljy = grid_parameters[11]

    # Dimension parameters
    nv, nu, nw = point_x.shape

    #--------------------------------------------------------------------------
    #
    # Calculate the mass fluxes
    #
    #--------------------------------------------------------------------------
    # Mass flux across each "i" face of elements.
    # Also compute total mass flow rate across each "i" line
    for j in range(0, nv-1):
        for i in range(0,nu):
            flux_i_mass[j,i] = 0.5 *(  ( ro_vel_x[j,i,0] + ro_vel_x[j+1, i,0]) * dlix[j,i] + (ro_vel_y[j,i,0] + ro_vel_y[j+1,i,0]) * dliy[j,i])
            flow[i,0] = flow[i,0] + flux_i_mass[j,i]

    # Mass flux across each "j" face of elements
    for j in range(1, nv):
        for i in range(0, nu - 1):
            flux_j_mass[j,i] = 0.5 * ( (ro_vel_x[j,i,0] + ro_vel_x[j,i+1,0]) * dljx[j,i]  +  (ro_vel_y[j,i] + ro_vel_y[j,i+1]) * dljy[j,i] )

    # Set the mass fluxes through j=1 and j=nv as zero -- these are solid surfaces
    # not necessary to resolve the velocity parallel to the surfaces
    for i in range(0, nu - 1):
        flux_j_mass[0,i] = 0.0
        flux_j_mass[nv-1,i] = 0.0

    #--------------------------------------------------------------------------
    #
    # Calculate the fluxes of momentum
    #
    #--------------------------------------------------------------------------
    # momentum in "i"
    for j in range(0, nv - 1):
        for i in range(0, nu):
            flux_i_xmom[j,i] = 0.5 * ( flux_i_mass[j,i] * (vel_x[j,i,0] + vel_x[j+1,i,0]) + (pressure[j,i,0] + pressure[j+1,i,0]) * dlix[j,i] )
            flux_i_ymom[j,i] = 0.5 * ( flux_i_mass[j,i] * (vel_y[j,i,0] + vel_y[j+1,i,0]) + (pressure[j,i,0] + pressure[j+1,i,0]) * dliy[j,i] )

    # momentum in "j"
    for j in range(0, nv):
        for i in range(0, nu - 1):
            flux_j_xmom[j,i] = 0.5 * ( flux_j_mass[j,i] * (vel_x[j,i,0] + vel_x[j,i+1,0]) + (pressure[j,i,0] + pressure[j,i+1,0]) * dljx[j,i] )
            flux_j_ymom[j,i] = 0.5 * ( flux_j_mass[j,i] * (vel_y[j,i,0] + vel_y[j,i+1,0]) + (pressure[j,i,0] + pressure[j,i+1,0]) * dljy[j,i] )

    #--------------------------------------------------------------------------
    #
    # Calculate the fluxes of enthalpy
    #
    #--------------------------------------------------------------------------
    # Enthalpy in "i"
    for j in range(0, nv - 1):
        for i in range(0, nu):
            flux_i_enthalpy[j,i] = 0.5 * flux_i_mass[j,i] * (enthalpy_stag[j,i] + enthalpy_stag[j+1,i])

    # Enthalpy in "j"
    for j in range(0, nv):
        for i in range(0, nu - 1):
            flux_j_enthalpy[j,i] = 0.5 * flux_j_mass[j,i] * (enthalpy_stag[j,i] + enthalpy_stag[j,i+1])

    # Save all this new data!
    fluxes[0] = flux_i_mass
    fluxes[1] = flux_j_mass
    fluxes[2] = flux_i_xmom
    fluxes[3] = flux_j_xmom
    fluxes[4] = flux_i_ymom
    fluxes[5] = flux_j_ymom
    fluxes[6] = flux_i_enthalpy
    fluxes[7] = flux_j_enthalpy
    fluxes[8] = flow

    return fluxes

"""
subroutine sums the fluxes for each element, calculates the changes in these
variable "prop" and distributes them to the four corners of the element
"""
def sum_fluxes(iflux, jflux, prop, prop_start, delprop, frkut, step, areas):

    #---------------------------------------------------------
    # FORMAT
    #
    # iflux, jflux = np.zeros((nv, nu))
    # prop = np.zeros((nv, nu, nw ))
    # NOTE: step comes from "time-stepping"
    #---------------------------------------------------------

    # Allocate "local" memory
    nv, nu, nw = prop.shape
    store = np.zeros((nv, nu))

    # Set store to be total flux in prop
    for j in range(0, nv - 1):
        for i in range(0, nu - 1):
            totalflux = iflux[j,i] - iflux[j,i+1] + jflux[j,i] - jflux[j+1,i]
            store[j,i] = frkut * step[j,i] * totalflux * areas[j,i]

    # Distribute the change equally to the four interior corners of each cell.
    # Each interior grid points receieve 1/4 of the change from the four
    # cells adjacent to it.
    for j in range(1, nv - 1):
        for i in range(1, nu - 1):
            add = 0.25 * (store[j,i] + store[j-1,i-1] + store[j-1,i] + store[j,i-1])
            prop[j,i,0] = prop_start[j,i,0] + frkut*add

    # add the changes to the upper and lower boundaries -- these receieve half the
    # change from each of the two cells adjacent to them...
    for i in range(1, nu - 1):

        # add to nodes with j = 0
        add = 0.5 * ( store[0,i] + store[0,i-1] )
        prop[0,i,0] = prop_start[0,i,0] + frkut*add

        # add to nodes with j = nv
        add = 0.5 * (store[nv-1,i] + store[nv-1,i])
        prop[nv-1,i,0] = prop_start[nv-1,i,0] + frkut * add

    # Now add on the extra changes to the inlet and outlet
    for j in range(1, nv - 1):

        # add to nodes with i = nu
        add = 0.5 * (store[j, nu - 1] + store[j-1, nu - 1] )
        prop[j, nu - 1, 0] = prop_start[j, nu - 1, 0 ] + frkut*add

        # add to ntoes with i = 0
        add = 0.5 * (store[j, 0] + store[j-1, 0])
        prop[j, 0, 0] = prop_start[j, 0, 0] + frkut * add

    # For the node with i = 0, j = 0
    add = store[0,0]
    prop[0,0,0] = prop_start[0,0,0] + frkut * add

    # For the node with i = 0, j = nv - 1
    add = store[nv-1,0]
    prop[nv-1,0,0] = prop_start[nv-1,0,0] + frkut * add

    # For the node with i = nu - 1, j = 0
    add = store[0, nu - 1]
    prop[0, nu - 1, 0] = prop_start[0, nu - 1, 0] + frkut * add

    # For the node with i = nu -1, j = nv - 1
    add = store[nv-1, nu-1]
    prop[nv-1,nu-1,0] = prop_start[nv-1,nu-1,0] + frkut * add

    # Now save the changes in the primary variable as delprop
    for j in range(0, nv-1):
        for i in range(0, nu-1):
            delprop[j,i] = store[j,i]

    return prop, delprop
