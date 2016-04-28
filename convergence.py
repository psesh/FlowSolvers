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
def check_convergence(primary_variables, old_primary_variables, reference_values, nstep):

    # Unpack the primary flow variables
    ro = primary_variables[0]
    ro_vel_x = primary_variables[1]
    ro_vel_y = primary_variables[2]
    ro_energy = primary_variables[3]

    # Unpack the old primary flow variables
    old_ro = old_primary_variables[0]
    old_ro_vel_x = old_primary_variables[1]
    old_ro_vel_y = old_primary_variables[2]
    old_ro_energy = old_primary_variables[3]

    # Unpack the "reference values"
    ref_ro = reference_values[0]
    ref_ro_vel_x = reference_values[1]
    ref_ro_vel_y = reference_values[2]
    ref_ro_energy = reference_values[3]


    # Declare some initial values
    del_ro_max = 0.0
    del_ro_vel_x_max = 0.0
    del_ro_vel_y_max = 0.0
    del_ro_energy_max = 0.0

    # Average values
    del_ro_avg = 0.0
    del_ro_vel_x_avg = 0.0
    del_ro_vel_y_avg = 0.0
    del_ro_energy_avg = 0.0


    imax = 0 # These are indices!
    jmax = 0

    #[imax, jmax] is the grid point where the change in rovx is a max
    for j in range(0, nv):
        for i in range(0, nu):

            # 1 ) Density
            delta = np.abs(ro[j,i,0] - old_ro[j, i, 0])
            if( delta > delromax ):
                del_ro_max = delta
                del_ro_avg = del_ro_avg + delta

            # 2 ) Density-X-velocity
            delta = np.abs(ro_vel_x[j,i,0] - old_ro_vel_x[j,i,0])
            if(delta > del_ro_vel_x_max):
                del_ro_vel_x_max = delta
                imax = i
                jmax = j
            del_ro_vel_x_avg = del_ro_vel_x_avg + delta

            # 3 ) Density-Y-velocity
            delta = np.abs(ro_vel_y[j,i,0] - old_ro_vel_y[j,i,0])
            if(delta > del_ro_vel_y_max):
                del_ro_vel_y_max = delta
            del_ro_vel_y_avg = del_ro_vel_y_avg + delta

            # 4 ) Density-Energy
            delta = np.abs(ro_energy[j,i,0] - old_ro_energy[j,i,0] )
            if(delta > del_ro_energy_max):
                del_ro_energy_max = delta
            del_ro_energy_avg = del_ro_energy_avg + delta


    # Compute the average changes
    del_ro_avg = del_ro_avg / ncells / reference_ro
    del_ro_vel_x_avg = del_ro_vel_x_avg / ncells / ref_ro_vel_x
    del_ro_vel_y_avg = del_ro_vel_y_avg / ncells / ref_ro_vel_y
    del_ro_energy_avg = del_ro_energy_avg / ncells /  ref_ro_energy
    del_ro_vel_x_max = del_ro_vel_x_max / ref_ro_vel_x
    del_ro_vel_y_max = del_ro_vel_y_max / ref_ro_vel_y

    emax = np.max([del_ro_vel_x_max, del_ro_vel_y_max])
    eavg = np.average([del_ro_vel_x_avg, del_ro_vel_y_avg])

    # Store the current values as the old ones
    for j in range(0, nv):
        for i in range(0, nu):
            old_ro[j,i,0] = ro[j,i,0]
            old_ro_vel_x[j,i,0] = old_ro_vel_x[j,i,0]
            old_ro_vel_y[j,i,0] = old_ro_vel_y[j,i,0]
            old_ro_energy[j,i,0] = old_ro_energy[j,i,0]

    # Pack these values!
    old_primary_variables[0] = old_ro
    old_primary_variables[1] = old_ro_vel_x
    old_primary_variables[2] = old_ro_vel_y
    old_primary_variables[3] = old_ro_energy

    # On-screen output
    print("------ Time Step Number ----- %i" %nstep)
    print("Emax = %d , at imax = %d , at jmax = %d , eavg = %d " %(emax, imax, jmax, eavg) )

    return old_primary_variables
