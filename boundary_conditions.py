#!/usr/bin/python
import numpy as np
import os


"""
Subroutine that calculates the secondary flow variables from the primary ones at
each grid point and then returns the computed secondary values


    Inputs:
        primary_variables
        secondary_variables

    Outputs:
        secondary_variables
"""
def set_other_variables(primary_variables, secondary_variables, boundary_conditions, grid_parameters):

    # Get the shape from the grid parameters
    point_x = grid_parameters[0]
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

    # "Unpack" the boundary conditions
    rgas = boundary_conditions[0]
    gamma = boundary_conditions[1]
    pressure_stag_inlet = boundary_conditions[2]
    temp_stag_inlet = boundary_conditions[3]
    alpha_1 = boundary_conditions[4]
    pressure_static_exit = boundary_conditions[5]

    # Additional flow properties
    cp = rgas * gamma / (gamma - 1.0)
    cv = cp / (gamma * 1.0)

    # Compute the new secondary values!
    for j in range(0, nv):
        for i in range(0, nu):
            vel_x[j,i,0] = ro_vel_x[j,i,0] /(1.0 * ro[j,i,0])
            vel_y[j,i,0] = ro_vel_y[j,i,0] /(1.0 * ro[j,i,0])
            enthalpy_stag[j,i,0] = cp * temp_stag_inlet
            kinetic_energy = 0.5 * (vel_x[j,i,0] * vel_x[j,i,0] + vel_y[j,i,0] * vel_y[j,i,0])
            temperature = (enthalpy_stag[j,i,0] - kinetic_energy ) / (1.0 * cp)
            pressure[j,i,0] = ro[j,i,0] * rgas * temperature


    # Packing the new values!
    secondary_variables[0] = vel_x
    secondary_variables[1] = vel_y
    secondary_variables[2] = pressure
    secondary_variables[3] = enthalpy_stag

    return secondary_variables

"""
Subroutine applies the boundary conditions that pressure = pressure_static_exit
at the exit. At the inlet boundary, the change in the density is relaxed.

    Inputs:
        nstep
        primary_variables
        secondary_variables
        boundary_conditions
        grid_parameters

    Outputs:
        primary_variables
        secondary_variables
"""
def apply_boundary_conditions(nstep, primary_variables, secondary_variables, boundary_conditions, grid_parameters):

    # Grid parameter for loops!
    point_x = grid_parameters[0]
    nv, nu, nw = point_x.shape

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

    # "Unpack" the boundary conditions
    rgas = boundary_conditions[0]
    gamma = boundary_conditions[1]
    pressure_stag_inlet = boundary_conditions[2]
    temp_stag_inlet = boundary_conditions[3]
    alpha_1 = boundary_conditions[4]
    pressure_static_exit = boundary_conditions[5]

    # allocate memory for "roinlet":
    ro_inlet = np.zeros((nv, 1))

    # As the inlet condition may become unstable, we relax the changes in
    # the inlet density by a factor of "rfin"
    density_fac = 0.25
    density_fac_one = 1 - density_fac
    ro_stag_inlet = pressure_stag_inlet / (rgas * temp_stag_inlet)
    cp = rgas * gamma / (gamma - 1.0)
    cv = cp / (gamma * 1.0)

    for j in range(0, nv):
        if(nstep == 1):
            ro_inlet[j,0] = ro[j,0,0]
        else:
            ro_inlet[j,0] = density_fac * ro[j,0,0] + density_fac_one * ro_inlet[j,0]

        if(ro_inlet[j,0] > 0.9999 * ro_stag_inlet):
            ro_inlet[j,0] = 0.9999 * ro_stag_inlet

        pressure[j,0,0] = pressure_stag_inlet * (ro_inlet[j,0] / ro_stag_inlet) ** gamma
        temperature = pressure[j,0,0] / (rgas * ro_inlet[j,0])
        velocity = np.sqrt(2.0 * cp * (temp_stag_inlet - temperature))
        vel_x[j,0,0] = velocity * np.cos(alpha_1)
        vel_y[j,0,0] = velocity * np.sin(alpha_1)
        ro_vel_x[j,0,0] = ro_inlet[j,0] * vel_x[j,0,0]
        ro_vel_y[j,0,0] = ro_inlet[j,0] * vel_y[j,0,0]
        ro_energy[j,0,0] = ro_inlet[j,0] * (cv * temperature + 0.5 * (velocity ** 2) )
        enthalpy_stag[j,0,0] = cp * temp_stag_inlet
        pressure[j, nu - 1, 0] = pressure_static_exit


    # Now that we are done, return the primary and secondary flow variables!
    primary_variables[0] = ro
    primary_variables[1] = ro_vel_x
    primary_variables[2] = ro_vel_y
    primary_variables[3] = ro_energy
    secondary_variables[0] = vel_x
    secondary_variables[1] = vel_y
    secondary_variables[2] = pressure
    secondary_variables[3] = enthalpy_stag

    return primary_variables, secondary_variables
