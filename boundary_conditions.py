#!/usr/bin/python
import numpy as np
import os

"""
Subroutine applies the boundary conditions that pressure = pressure_static_exit
at the exit. At the inlet boundary, the change in the density is relaxed.
"""
def apply_boundary_conditions(primary_variables, secondary_variables, boundary_conditions, grid_parameters):

    # Grid parameter for loops!
    point_x = grid_parameters[0]
    nv, nu, nw = point_x.shape

    # "Unpack" the primary flow variables
    ro_vel_x = primary_variables[1]
    ro_vel_y = primary_variables[2]

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
    pressure_static_exit = boundary_conditions[5]

    # As the inlet condition may become unstable, we relax the changes in
    # the inlet density by a factor of "rfin"
    density_fac = 0.25
    density_fac_one = 1 - density_fac
    ro_stag_inlet = pressure_stag_inlet / (rgas * temp_stag_inlet)

    for j in range(0, nv)








    rfin = 0.25
    rfin1 = 1.0 - rfin
    ro_stag_inlet = pressure_stag_inlet / (rgas * )
    pressure_stag_inlet =
    rostagin = pstagin / (rgas * tstagin)

    for i in range(0, nj-1):
        if(nstep == 1):
            roinlet[j] = ro[0,j]
        else:
            roinlet[j] = rfin * ro[0,j] + rfin1 * roinlet[j]

    if(roinlet[j] > 0.999 * rostagin):
        roinlet[j] = 0.999 * rostagin

    """ Insert code here to calculate P[1,j], rovx[1,j], rovy[1,j]
    and roe[1,j] from roinlet[j], pstagin, tstagin and alpha1. also set
    vx[1,j], vy[1,j] and hstag[1,j]
    """

    p[0,j] = pstagin * (roinlet[j] / rostagin) ** gamma
    t = P[0,j] / (rgas * roinlet[j])
    vel = (2 * cp * (tstagin - t)) ** 0.5

    vx[1,j] = vel * np.cos(alpha1)
    vy[1,j] = vel * np.sin(alpha1)

    rovx[1,j] = roinlet[j] * vx[1,j]
    rovy[1,j] = roinlet[j] * vy[1,j]

    roe[1,j] = roinlet[j] * vy[1,j]
    hstag[1,j] = cp * tstagin
    p[ni,j] = pdown
