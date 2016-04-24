#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
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
def smooth(property, corr_prop, boundary_conditions, nu, nv):

    # constants
    smooth_fac_input = boundary_conditions[7]
    smooth_fac_minus_one = 1.0 - smooth_fac_input
    fcorrection = 0.95

    for i in range(0, nu):
        ip1 = i + 1
        if i == nu - 1 :
            ip1 = nu - 1
        if i == 1 :
            im1 = 1

        for j in range(2, )
