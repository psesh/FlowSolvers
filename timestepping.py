#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy as np



def set_timestep(primary_variables, secondary_variables, boundary_conditions):

    # Subroutine sets the length of the timestep based on the stagnation speed of sound
    # and the minimum lenght scale of any element. The timestep must be called "delta_t"
    # An assumption that the maximum flow speed will equal to the "Astag" is also made.
    # This will be pessismistic for subsonic flows but may be optimic for supersonic flows.
    # cfl number was an

    # Unpack the primary f& secondary flow variables
    vel_x = secondary_variables[0]
    vel_y = secondary_variables[1]
    ro = primary_variables[0]
    pressure = secondary_variables[2]

    gamma = boundary_conditions[1]
    cfl = boundary_conditions[6]

    boundary_conditions =
    ## ---> [rgas, gamma, pressure_stag_inlet, temp_stag_inlet, alpha_1, pressure_static_exit, cfl ,  smooth_fac_input, nsteps, conlim_in ]


    # Velocity magnitude
    velocity_magnitude = np.zeros((nu, nv))

    # Setting the velocity magnitude!
    for j in range(0, nv):
        for i in range(0, nu):
            velocity_magnitude[j,i] = np.sqrt( vel_x[j,i,0] * vel_x[j,i,0] + vel_y[j,i,0] * vel_y[j,i,0] )

    for j in range(0, nv - 1):
        for i in range(0, nu - 1):
            A_1 = np.sqrt(gamma * pressure[j,i,0]/ro[j,i,0]
            A_2 = 
