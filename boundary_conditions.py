#!/usr/bin/python
import numpy as np
import os

"""
This routine applies the boundary conditions that P = Pdown at i = ni. At the
inlet boundary the change in density is relaxed to obtain roinlet(j) which is then
used to obtain other properties at the inlet assuming isentropic flow from
stagnation conditions, "po and to" together with the specified inlet flow angle
alpha1. Because the inlet condition may become unstable, it is safer to relax
the changes in inlet density by a factor "rfin = 0.25" as set below. Reduce this
further if the inlet becomes unstable. Also worth checking roinlet.
"""
def apply_boundary_conditions():

    rfin = 0.25
    rfin1 = 1.0 - rfin
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
    
    