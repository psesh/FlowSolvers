#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import numpy as np


"""
    Subroutine sets the length of the timestep based on the stagnation speed of sound
    and the minimum lenght scale of any element. The timestep must be called "delta_t"
    An assumption that the maximum flow speed will equal to the "Astag" is also made.
    This will be pessismistic for subsonic flows but may be optimic for supersonic flows.
    cfl number was an
"""
def set_timestep(primary_variables, secondary_variables, boundary_conditions, grid_parameters):

    # Unpack the primary and secondary flow variables
    point_x = grid_parameters[0]
    vel_x = secondary_variables[0]
    vel_y = secondary_variables[1]
    ro = primary_variables[0]
    pressure = secondary_variables[2]
    gamma = boundary_conditions[1]
    cfl = boundary_conditions[6]
    dmin = grid_parameters[12]
    nv, nu, nw = point_x.shape

    # Velocity magnitude
    velocity_magnitude = np.zeros((nv, nu))

    # Time step parameter
    step = np.zeros((nv, nu))

    # Setting the velocity magnitude!
    for j in range(0, nv):
        for i in range(0, nu):
            velocity_magnitude[j,i] = np.sqrt( vel_x[j,i,0] * vel_x[j,i,0] + vel_y[j,i,0] * vel_y[j,i,0] )

    for j in range(0, nv - 1):
        for i in range(0, nu - 1):
            a1 = np.sqrt(gamma * pressure[j,i,0]/ro[j,i,0] )
            a2 = np.sqrt(gamma * pressure[j,i+1,0]/ro[j,i+1,0] )
            a3 = np.sqrt(gamma * pressure[j+1,i, 0]/ro[j,i+1,0] )
            a4 = np.sqrt(gamma * pressure[j+1,i+1,0]/ro[j+1,i+1,0])
            aaverage = 0.25 * (a1 + a2 + a3 + a4)

            velocity1 = velocity_magnitude[j,i]
            velocity2 = velocity_magnitude[j,i+1]
            velocity3 = velocity_magnitude[j+1,i]
            velocity4 = velocity_magnitude[j+1,i+1]
            velaverage = 0.25 * (velocity1 + velocity2 + velocity3 + velocity4)

            step[j,i] = cfl * dmin[j,i]/(aaverage + velaverage)

    # Manually place values for step[j,i] for the last row and column
    for i in range(0, nu - 1):
        aaverage = np.sqrt(gamma * pressure[nv - 1 , i] / ro[nv - 1 , i])
        velaverage = velocity_magnitude[nv - 1, i]
        step[nv - 1, i] = cfl * dmin[nv - 1, j] / (aaverage + velaverage)

    for j in range(0, nv - 1):
        aaverage = np.sqrt(gamma * pressure[j, nu - 1] / ro[j, nu - 1])
        velaverage = velocity_magnitude[j, nu - 1]
        step[j, nu - 1] = cfl * dmin[j, nu - 1] / (aaverage + velaverage)

    aaverage = np.sqrt(gamma * pressure[nv - 1, nu - 1] / ro[nv - 1, nu - 1])
    step[nv - 1, nu - 1] = cfl * dmin[nv - 1, nu - 1]/ (velocity_magnitude[nv - 1, nu - 1] + aaverage)

    return step
