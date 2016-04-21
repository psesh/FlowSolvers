#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
from evtk.hl import gridToVTK, pointsToVTK
from geometry import create_grid
import numpy as np

# Generate an initial flow solution based off the flow conditions and the geometry
def initial_guess():
    
    nu, nv, nw = 120, 60, 1 
    #----------------------------------------------
    #
    # Here we initialize the flow variables!
    #
    #-----------------------------------------------    
    ro = np.zeros((nv, nu, nw)) # Density
    ro_vx = np.zeros((nv, nu, nw)) # Density * Velocity in "x" direction
    ro_vy = np.zeros((nv, nu, nw)) # Density * Velocity in "y" direction
    pressure = np.zeros((nv, nu, nw)) # static pressure
    enthalpy_stagnation = np.zeros((nv, nu, nw)) # stagnation enthalpy
    ro_energy = np.zeros((nv, nu, nw)) # Density * energy
    
    # Fluid / gas constants
    gamma = 1.4 # gas specific heat ratio
    rgas  = 287.5 # gas constant J/kg K
    
    # Dual-for loops for initialization!
    for k in range(0, nw):
        for j in range(0, nv):
            for i in range(0, nu):
                ro[j,i,k] = 1.2
                ro_vx[j,i,k] = 100.0 * j/(nv * 1.0)
                ro_vy[j,i,k] = 0.0
                pressure[j,i,k] = 100000 * (0.9 + 0.1 * j/(nv * 1.0))
                enthalpy_stagnation[j,i,k] = 300000.0
                ro_energy[j,i,k] = pressure[j,i,k] / (gamma - 1.0)
    
    # Calculate the reference values for checking convergence
    ncells = nv * nu * nw
    jmid = (1 + nv) / (2.0)
    #ro_in = pressure_stagnation / rgas / 
    
    #pressure = np.random.rand(ncells).reshape((nv-1, nu-1, nw-1))
    #temp = np.random.rand(npoints).reshape((nv, nu, nw))

    # Converting from points to VTK readable format!
    point_x, point_y, point_z, areas = create_grid(nu, nv, nw)
    gridToVTK("./structured", point_x, point_y, point_z, pointData={"pressure": pressure, "density": ro })

initial_guess()