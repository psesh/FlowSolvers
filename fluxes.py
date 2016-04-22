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
    ro = primary_variables[0]
    ro_vel_x = primary_variables[1] 
    ro_vel_y = primary_variables[2]  
    ro_energy = primary_variables[3] 
    
    # "Unpack" the secondary flow variables
    vel_x = secondary_variables[0]
    vel_y = secondary_variables[1]
    pressure = secondary_variables[2] 
    enthalpy_stag = secondary_variables[3] 
    
    # "Unpack" the grid parameters
    point_x = grid_parameters[0] 
    point_y = grid_parameters[1] 
    point_z = grid_parameters[2] 
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
    fluxes[0] = flux_i_mass, fluxes[1] = flux_j_mass, fluxes[2] = flux_i_xmom
    fluxes[3] = flux_j_xmom, fluxes[4] = flux_i_ymom, fluxes[5] = flux_j_ymom
    fluxes[6] = flux_i_enthalpy, fluxes[7] = flux_j_enthalpy
    fluxes[8] = flow
    
    return fluxes
"""
subroutine sums the fluxes for each element, calculates the changes in these
variable "prop" and distributes them to the four corners of the element
"""
def sum_fluxes(iflux, jflux, prop, delprop, frkut):
    
    #---------------------------------------------------------
    # FORMAT
    #
    # iflux, jflux = np.zeros((nv, nu))
    # prop = np.zeros((nv, nu, nw ))
    # 
    #---------------------------------------------------------
    
    # Allocate "local" memory
    store = np.zeros((nv, nu))
    
    for j in range(0, nv - 1):
        for i in range(0, nu - 1)"
            totalflux = 
    
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
