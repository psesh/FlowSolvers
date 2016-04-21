#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
from evtk.hl import gridToVTK, pointsToVTK
from geometry import create_grid
import numpy as np

# Generate an initial flow solution based off the flow conditions and the geometry
def initial_guess():
    
    # Geometry parameters:
    nu, nv, nw = 120, 60, 1
    
    
    #--------------------------------------------------------------
    #
    # Read in the boundary conditions file and save variables
    #
    #--------------------------------------------------------------
    
    # Read in boundary condition file!
    bcfile = open('boundary_conditions.txt', 'r')
    input_data = bcfile.readlines()
    bcfile.close()
    
    # Skip the header and proceed with the initialization!
    line_2 = input_data[1].split()
    rgas = np.float(line_2[2])
    
    line_3 = input_data[2].split()
    gamma = np.float(line_3[2])

    line_4 = input_data[3].split()
    pressure_stag_inlet = np.float(line_4[2])

    line_5 = input_data[4].split()
    temp_stag_inlet = np.float(line_5[2])

    line_6 = input_data[5].split()
    alpha_1 = np.float(line_6[2])

    line_7 = input_data[6].split()
    pressure_static_exit = np.float(line_7[2])
    
    line_8 = input_data[7].split()
    cfl = np.float(line_8[2])
    
    line_9 = input_data[8].split()
    smooth_fac_input = np.float(line_9[2])
    
    line_10 = input_data[9].split()
    nsteps = np.float(line_10[2])
    
    line_11 = input_data[10].split()
    conlim_in = np.float(line_11[2])

    #----------------------------------------------
    #
    # Setup other variables!
    #
    #----------------------------------------------
    emax = 1000000
    eavg = emax
    cp = rgas * gamma / (gamma - 1.0)
    cv = cp / (gamma * 1.0)
    gamma_factor = (gamma - 1.0) / (gamma * 1.0)
    smooth_fac = smooth_fac_input * cfl
    conlim = conlim_in * cfl
    alpha_1 = alpha_1 * np.pi/(180.0)
    
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
    

    # Dual-for loops for initialization!
    for k in range(0, nw):
        for j in range(0, nv):
            for i in range(0, nu):
                ro[j,i,k] = 1.2
                ro_vx[j,i,k] = 100.0 * j/(nv * 1.0)
                ro_vy[j,i,k] = 0.0
                pressure[j,i,k] = 100000 * (0.9 + 0.1 * i/(nu * 1.0))
                enthalpy_stagnation[j,i,k] = 300000.0
                ro_energy[j,i,k] = pressure[j,i,k] / (gamma - 1.0)
    
    # Calculate the reference values for checking convergence
    ncells = nv * nu * nw
    jmid = (1 + nu) / (2.0)
    ro_in = pressure_stag_inlet / rgas / temp_stag_inlet
    ref_ro = (pressure_stag_inlet - pressure_static_exit) / rgas / temp_stag_inlet
    ref_temp = temp_stag_inlet * (pressure_static_exit / pressure_stag_inlet) ** gamma_factor
    ref_velocity = np.sqrt(2 * cp * (temp_stag_inlet - ref_temp))
    ref_ro_vx = ro_in * ref_velocity
    ref_ro_vy = ref_ro_vx
    ref_ro_energy = ro_in * cv * (temp_stag_inlet - ref_temp)
    
    #pressure = np.random.rand(ncells).reshape((nv-1, nu-1, nw-1))
    #temp = np.random.rand(npoints).reshape((nv, nu, nw))

    # Converting from points to VTK readable format!
    point_x, point_y, point_z, xlow, ylow, xhigh, yhigh, areas = create_grid(nu, nv, nw)
    gridToVTK("./structured", point_x, point_y, point_z, pointData={"pressure": pressure, "density": ro })

# More appropriate estimate of flow conditions!
def refined_flow_estimate():

    # Subroutine where we make an initial guess of the primary flow variables:
    # ro, ro_vx, ro_vy, ro_energy
    # The guess does not need to be very accuracy, but the better it is the 
    # faster the program will converge. Values to the primary flow variables
    # will be given at each grid point in this subroutine

    # Allocate memory
    aflow = np.zeros((nu, 1)) # length of each "i" between points
    ro_guess = np.zeros((nu, 1))
    vel_guess = np.zeros((nu, 1))
    
    # Set maximum values
    mach_limit = 1.0
    temp_limit = temp_stag_inlet / (1.0 + 0.5 * (gamma - 1.0) * mach_limit * mach_limit)
    
    # Work out the length of each "i" line between grid points, "i, 1" and "i, nj"
    # and call it AFLOW[I]
    for i in range(0, nu):
        aflow[i,0] = np.sqrt( (xhigh[i,0] - xlow[i,0])**2 + (yhigh[i,0] - ylow[i,0])**2 )
        
    # Make an initial guess of the density and the velocity at the exit by assuming 
    # isentropic flow conditions. 
    ro_stag_in = pressure_stag_inlet / (rgas * temp_stag_inlet)
    ro_exit = ro_stag_in * (pressure_static_exit / pressure_stag_inlet) ** (1.0/gamma)
    temp_static_exit = pressure_static_exit / (rgas * ro_exit)
    vel_exit = np.sqrt(2 * cp * (temp_stag_inlet - temp_static_exit))
    mflow = ro_exit * aflow[nu,0] * vel_exit
    
    # Now estimate the velocity and density at every "i" line. Call the velocity
    # V_guess[i] and the density ro_guess[i]. First assume that the density is constant
    # and that the flow is perpendicular to the "i"=constant lines and hence occupies
    # the area aflow[i].
    # Hence use continuity to estimate the flow velocity v_guess[i]
    # Use this velocity to calculate the static temperature assuming that the 
    # stagnation temperature is constant
    # Check that this temperature is not less than TLIM and set = TLIM if it is
    # Next use this temperature and isentropic flow to obtain a better estimate
    # of this density, ro_guess[i]
    # Use this density and continuity to obtain a better estimate of the velocity
    # and set = V_guess[i]
    
    for i in range(0, nu-1):
        ro_guess[i,0] = ro_exit
        vel_guess[i,0] = mflow / (ro_guess[i,0] * aflow[i,0])
        temp_static = temp_stag_inlet - 0.5 *  (vel_guess[i,0]**2)/(cp)
        
        # Check if temp_static is within limits
        if temp_static < temp_limit:
            temp_static = temp_limit
        
        ro_guess[i,0] = 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
        
initial_guess()