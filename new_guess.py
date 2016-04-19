# -*- coding: utf-8 -*-
import numpy as np
"""
Created on Tue Apr 19 08:24:08 2016

@author: pseshadr
"""

def new_guess():
    
    jmid = nj / 2
    tdown = tstagin * (pdown / pstagin) ** fga
    vdown = np.sqrt( 2 * cp * (tstagin - tdown ))
    rodown = pdown / rgas / tdown
    pinlet = 55000
    tinlet = tstagin * (pinlet / pstagin) ** fga
    vinlet = np.sqrt(2 * cp * (tstagin - tinlet ))
    roin = pinlet / rgas / tinlet
    
    for j in range(0, nj):
        for i in range(0, ni - 2):
            dx = x[i + 1, jmid] - x[i, jmid]
            dy = y[i + 1, jmid] - y[i, jmid]
            ds = np.sqrt(dx * dx + dy * dy)
            vlocal = vinlet + (vdown - vinlet) * (i - 1.0)/(ni - 1.0)
            rolocal = roin + (rodown - roin) * (i - 1.0)/(ni - 1.0)
            tlocal = tinlet + (tdown - tinlet) * (i = 1.0)/(ni - 1.0)
            
            xvel = vlocal * dx / ds
            yvel = vlocal * dy / ds
            rovx[i,j] = rolocal * xvel
            rovy[i,j] = rolocal * yvel
            ro[i,j] = rolocal
            roe[i,j] = rolocal * (cv * tlocal  + 0.5 * vlocal * vlocal)
            
        rovx[ni, j] = rovx[ni - 1, j]
        rovy[ni, j] = rovy[ni - 1, j]
        ro[ni, j] = ro[ni - 1, j]
        roe[ni, j] = roe[ni - 1, j]
        
        