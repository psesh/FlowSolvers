#!/usr/bin/python
import numpy as np

"""
Subroutine that checks the changes in all the primary variables over the last 5
steps and writes a short output summary on screen and to a file!

    Inputs:
        primary flow variables

    Outputs:
        ??
"""
def check_convergence():

    delromax = 0.0
    delrovxmax = 0.0
    delrovymax = 0.0
    delroemax = 0.0
    delroavg = 0.0
    delrovxavg = 0.0
    delrovyavg = 0.0
    delroeavg = 0.0
    imax = 0.0
    jmax = 0.0

    #[imax, jmax] is the grid point where the change in rovx is a max
    for i in range(0, ni):
        for j in range(0, nj):

            delta = np.abs(ro[i,j] - ro_old[i,j] )
            if( delta > delromax ):
                delromax = delta

            delta = np.abs(rovx[i,j] - rovx_old[i,j])
            if(delta > delrovxmax ):
                imax = i
                jmax = j

            diffrovx[i,j] = delta
            delrovxavg = delrovxavg + delta
            delta = np.abs(rovy[i,j] - rovy_old[i,j])
            if( delta > delrovymax):
                delrovymax = delta
            delrovyavg = delrovyavg + delta

            delta = np.abs(roe[i,j] - roe_old[i,j])
            if( delta > delroemax):
                delroemax = delta
            delroeavg = delroeavg + delta

    # Compute the average changes
    delroavg = delroavg / ncells / ref_ro
    delrovxavg = delrovxavg / ncells / ref_rovx
    delrovyavg = delrovyavg / ncells / ref_rovy
    delroeavg = delroeavg / ncells / ref_roe
    delrovxmax = delrovxmax / ref_rovx
    delrovymax = delrovymax / ref_rovy

    emax = amax1(delrovymax, delrovymax)
    eavg = amax1(delrovxavg, delrovyavg)

    # store the maximum change in rovx and emax to be printed out
    # save the current values of the primary variables as prop_old values
    # for use in the next convergence check
    for i in range(0, ni - 1):
        for j in range(0, nj - 1):
            ro_old[i,j] = ro[i,j]
            rovx_old[i,j] = rovx[i,j]
            rovy_old[i,j] = rovy[i,j]
            roe_old[i,j] = roe[i,j]

    # Write the average changes in the primary variables to unit 3
    # for use in the convergence plotting program
    write_to_file(delroavg, delrovxavg, delrovyavg, delroeavg)

    # write a short output summary which will be written on screen
    flow_ratio = flow[ni - 1] / flow[0]
    print("------ Time Step Number ----- %i" %nstep)
    print("Emax = %d , at imax = %d , at jmax = %d , eavg = %d " %(emax, imax, jmax, eavg) )
    print("Inlet flow = %d , outlet to inlet flow ratio = %d ", %(flow[0], flow_ratio))
