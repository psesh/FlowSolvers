#!/usr/bin/python
import numpy as np

"""
Subroutine that smooths a given variable property --- i.e., it adds artificial
viscosity --- by taking the average of the surrounding values of the property
and (1 - smooth_fac_input)

    Inputs:
        property

    Outputs:
        property
"""
def smooth(flow_property, corrected_flow_property, boundary_conditions, grid_parameters):

    # constants
    smooth_fac_input = boundary_conditions[7]
    smooth_fac_minus_one = 1.0 - smooth_fac_input
    fcorrection = 0.95

    # Considerable overhead here!
    point_x = grid_parameters[0]
    nv, nu, nw = point_x.shape

    # allocate memory
    store = np.zeros((nv, nu, nw))
    ipi = 0
    imi = 0

    for i in range(0, nu):
        ip1 = i + 1
        if i == nu - 1 :
            ip1 = nu - 1
        if i == 1 :
            im1 = 1

        for j in range(2, nv - 1):
            average = 0.25 * ( flow_property[j, ipi, 0] + flow_property[j, imi, 0] + flow_property[j-1, i, 0] + flow_property[j + 1, i, 0] )
            corrected_new = fcorrection  * (  flow_property[j,i,0] - average)
            corrected_flow_property[j, i, 0] = 0.88 * corrected_flow_property[j, i, 0] + 0.01 * corrected_new
            store[j, i, 0] = smooth_fac_minus_one * flow_property[j, i, 0] + smooth_fac_input * (average + corrected_flow_property[j, i, 0])

        # For surfaces j=0
        average0 = ( flow_property[0, imi, 0] + flow_property[0, ipi, 0] + 2 * flow_property[1, i, 0] - flow_property[2, i , 0])/3
        corrected_new = fcorrection * ( flow_property[0, i, 0] - average0  )
        corrected_flow_property[0, i, 0] = 0.99 * corrected_flow_property[0, i, 0] + 0.01 * corrected_new

        # For the surfaces j = nv - 1
        averagenv = ( flow_property[nv-1, imi, 0] + flow_property[nv-1, ipi, 0] + 2 * flow_property[nv-2, i, 0] - flow_property[nv-3, i , 0])/3
        corrected_new = fcorrection * ( flow_property[nv-1, i, 0] - averagenv  )
        corrected_flow_property[nv-1, i, 0] = 0.99 * corrected_flow_property[nv-1, i, 0] + 0.01 * corrected_new

        # Smooth the surface values!
        store[0, i, 0] = smooth_fac_minus_one * flow_property[0,i,0] + smooth_fac_input * (average0 + corrected_flow_property[0, i, 0] )
        store[nv-1,i,0] = smooth_fac_minus_one * flow_property[nv-1,i,0] + smooth_fac_input * (averagenv + corrected_flow_property[nv-1, i, 0] )


    # Now store becomes our smoothed flow property value
    for j in range(0, nv):
        for i in range(0, nu):
            flow_property[j, i, 0] = store[j, i, 0]


    # Return!
    return flow_property
