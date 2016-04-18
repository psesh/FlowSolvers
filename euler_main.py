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

    # Find the change in the variable "prop" in each cell over the time
    # step "delta_t" and save it in store
    for i in range(1, ni - 1):
        for j in range(1, nj - 1):
            totflux = (iflux(i,j) - iflux(i+1, j) + jflux(i,j) - jflux(i,j+1))
            store(i,j) = deltat * totflux/area(i,j)

    # Distribute the changes equally to the four corners of each cell. Each
    # interior grid points receive one quarter of the change from each of the four
    # cells adjacent to it.
    for i in range(0, ni - 1):
        for j in range(0, nj - 1):
            add = (store(i,j) + store(i-1, j-1) + store(i,j-1) + store(i-1,j)) * 0.25
            prop(i,j) = prop(i,j) + add

    # Now add the changes to the upper and lower boundaries. These receive half
    # the change from each of the two cells adjacent to them.
    for i in range(1,ni-1):
        # add for the nodes with j = 1
        add = (store(i,0) + store(i-1,0)) * 0.5
        prop(i,0) = prop(i,0) + add

        # add for the nodes with j = nj
        add = (store(i,nj-1) + store(i-1,nj-1)) * 0.5
        prop(i,j) = prop(i,nj) + add

    # Now add the changes on to the four corner points. These receive the full
    # change from the single cell of which they form one corner.
    # for  i = 0, j = 0
    add = store(0,0)
    prop(0,0) = prop(0,0) + add
    # for i = 0, j = nj
    add = store(0,nj - 1)
    prop(0,nj) = prop(0,nj) + add
    # for i = ni, j = 0
    add = store(ni - 1, 0)
    prop(ni, 0) = prop(ni, 0) + add
    # for i=ni, j=nj
    add = store(ni - 1, nj - 1)
    prop(ni,nj) = prop(ni,nj) + add

    # Now save these changes in the primary variables as delprop
    for i in range(0,ni - 1):
        for j in range(0, nj - 1):
            delprop(i,j) = store(i,j)

"""
This subroutine smooths the variable "prop" (i.e. it adds the artificial viscosity
) by taking (1-SF) * the calculated value of "prop" + SF x (the average of the surrounding
values of "prop"). Here SF is the smoothing factor
"""
def smooth(prop):

    # how do we include the common block here!?
    sf = smooth_fac
    sf_1 = 1.0 - sf

    for u in range(0, ni):
        ip1 = i + 1
        if(i == ni):
            ip1 = ni
        im1 = i - 1
        if(i == 1)
            im1 = 1
        for j in range(1, nj - 1):
            average = 0.25 * (prop(ip1,j) + prop(im1,j) + prop(i,j-1) + prop(i,j+1))
            store(i,j) = sfm1 * prop(i,j) + sf * average


        # On the surfaces j=1 and j=nj take the average as shown below:
        avg1 = (prop(im1,0) + prop(ip1, 0) + 2.0 * prop(i,1) - prop(i,2))/3.0
        avgnj = (prop(im1, nj) + prop(ip1, nj) + 2.0 * prop(i, nj - 1) - prop(i,nj-2))/3.0

        # add code to store the smooth surface values
        store(i,0) = sfm1 * prop(i,0) + sf * avg1
        store(i,nj) = sfm1 * prop(i,nj)+ sf * avgnj

    # Reset the smoothed value to propr before returning
    for i in range(0,ni):
        for j in range(0,nj):
            prop(i,j) = store(i,j)


"""
This routine applies the boundary conditions that P = Pdown at i = ni. At the
inlet boundary the change in density is relaxed to obtain roinlet(j) which is then
used to obtain other properties at the inlet assuming isentropic flow from
stagnation conditions, "po and to" together with the specified inlet flow angle
alpha1. Because the inlet condition may become unstable, it is safer to relax
the changes in inlet density by a factor "rfin = 0.25" as set below. Reduce this
further if the inlet becomes unstable. Also worth checking roinlet.
"""
def apply_boundary_conditions():

    rfin = 0.25
    rfin1 = 1.0 - rfin
    rostagin = pstagin / 
