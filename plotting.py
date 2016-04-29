#!/usr/bin/python
from evtk.hl import gridToVTK
import numpy as np

def plot_to_grid(grid_parameters, primary_variables, secondary_variables, number):


    # Get the shape from the grid parameters
    point_x = grid_parameters[0]
    point_y = grid_parameters[1]
    point_z = grid_parameters[2]
    nv, nu, nw = point_x.shape

    # "Unpack" the primary flow variables
    ro = primary_variables[0]
    ro_vel_x = primary_variables[1]
    ro_vel_y = primary_variables[2]
    ro_energy = primary_variables[3]

    # "Unpack" the existing secondary flow variables
    vel_x = secondary_variables[0]
    vel_y = secondary_variables[1]
    pressure = secondary_variables[2]
    enthalpy_stag = secondary_variables[3]

    # mach number
    velocity_magnitude = np.zeros((nv, nu, nw))
    for j in range(0, nv):
        for i in range(0, nu):
            velocity_magnitude[j,i,0] = np.sqrt( vel_x[j,i,0] * vel_x[j,i,0] + vel_y[j,i,0] * vel_y[j,i,0] )

    filename = "./output"+str(number)
    gridToVTK(filename, point_x, point_y, point_z, pointData={"pressure": pressure, "density": ro, "vel-x": vel_x, "vel-y": vel_y , "vel-mag": velocity_magnitude})
