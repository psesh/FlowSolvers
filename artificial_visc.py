# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 07:52:16 2016

@author: pseshadr
"""

def smooth(prop):

    # how do we include the common block here!?
    sf = smooth_fac
    sf_1 = 1.0 - sf

    for u in range(0, ni):
        ip1 = i + 1
        if(i == ni):
            ip1 = ni
        im1 = i - 1
        if(i == 1)
            im1 = 1
        for j in range(1, nj - 1):
            average = 0.25 * (prop(ip1,j) + prop(im1,j) + prop(i,j-1) + prop(i,j+1))
            store(i,j) = sfm1 * prop(i,j) + sf * average


        # On the surfaces j=1 and j=nj take the average as shown below:
        avg1 = (prop(im1,0) + prop(ip1, 0) + 2.0 * prop(i,1) - prop(i,2))/3.0
        avgnj = (prop(im1, nj) + prop(ip1, nj) + 2.0 * prop(i, nj - 1) - prop(i,nj-2))/3.0

        # add code to store the smooth surface values
        store(i,0) = sfm1 * prop(i,0) + sf * avg1
        store(i,nj) = sfm1 * prop(i,nj)+ sf * avgnj

    # Reset the smoothed value to propr before returning
    for i in range(0,ni):
        for j in range(0,nj):
            prop(i,j) = store(i,j)
