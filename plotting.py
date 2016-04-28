#!/usr/bin/python
from evtk.hl import gridToVTK


def plot_to_grid(grid_parameters, primary_variables, secondary_variables):


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

    gridToVTK("./output", point_x, point_y, point_z, pointData={"pressure": pressure, "density": ro, "velx": vel_x, "vely": vel_y })
