#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
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
def set_fluxes(primary_variables, secondary_variables, fluxes, boundary_conditions):
    
        
    # "Unpack" the fluxes
    flux_i_mass = fluxes[0], 
    flux_j_mass = fluxes[1], 
    flux_i_xmom = fluxes[2]
    flux_j_xmom = fluxes[3] 
    flux_i_ymom = fluxes[4]
    flux_j_ymom = fluxes[5]
    flux_i_enthalpy = fluxes[6]
    flux_j_enthalpy = fluxes[7]
    
    # "Unpack" the primary flow variables
        
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
            totflux = (iflux[i,j
            \) - iflux(i+1, j) + jflux(i,j) - jflux(i,j+1))
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
