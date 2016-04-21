#!/usr/bin/python
#from evtk.hl import gridToVTK, pointsToVTK
import numpy as np
import random as rnd

"""
    Objectives:
    1. Generate a 2D converging-diverging structured grid with prescribed equations for the inlet and exit
    2. Export the grid to paraview for visualization

    Pranay Seshadri
    April 19th, 2016
    ps583@cam.ac.uk

"""
# Compute the areas of grid -- we do this via the cross product formula!
def compute_areas(x, y):

    # We think of each cell as a quadrilateral. Our objective is to compute
    # its area by taking the cross product of AC with BD. We begin by defining
    # the vectors AC and BD for each cell and then computing its area!
    nv, nu, nw = x.shape

    areas = np.zeros((nv, nu))
    for j in range(0, nv - 1):
        for i in range(0, nu - 1):
            A = np.array((x[j,i,0] ,  y[j,i,0] ))
            B = np.array(( x[j,i+1,0], y[j,i,0] ))
            C = np.array(( x[j+1,i+1,0], y[j+1,i+1,0]) )
            D = np.array(( x[j+1,i,0], y[j+1,i,0] ) )
            AC = np.array((C[0] - A[0], C[1] - A[1], 0 ))
            BD = np.array((D[0] - B[0], D[1] - B[1], 0 ))
            areas[j,i] = 1/2 * np.linalg.norm(np.cross(AC, BD))
    
    return areas

# Create a grid!
def create_grid(nu, nv, nw):

    real_width = 10.0
    real_height = 2.0
    xc = real_width/2.0
    yc = -20.0
    r = np.sqrt(xc*xc + yc*yc)


    ncells = (nv - 1) * (nu - 1) * (nw - 1)
    npoints = nv * nu * nw

    # Storing data for each point in the grid!
    point_x = np.zeros((nv , nu , nw ))
    point_y = np.zeros((nv , nu , nw ))
    point_z = np.zeros((nv , nu , nw ))
    
    # Storing upper and lower boundary values only!
    xlow = np.zeros((nu, 1))
    ylow = np.zeros((nu, 1))
    xhigh = np.zeros((nu, 1))
    yhigh = np.zeros((nu, 1))
    
    # Inlet
    for i in range(0, nv):
        point_x[i,0,0] = 0.0
        point_y[i,0,0] = i * real_height/(nv - 1)
        point_z[i,0,0] = 0

    # Exit
    for i in range(0, nv):
        point_x[i,nu-1,0] = real_width
        point_y[i,nu-1,0] = i * real_height/(nv - 1.0)
        point_z[i,nu-1,0] = 0

    # Parameters to define the top and bottom walls
    theta1 = np.arccos(real_width/(2.0 * r))
    theta2 = np.arccos(-real_width/(2.0 * r))
    x1 = r * np.cos(theta1)
    y1 = r * np.sin(theta1)

    # Bottom
    for i in range(0, nu):
        point_x[0,i,0] = x1 - r*np.cos(theta1 + (theta2 - theta1)/(nu-1)*i)
        point_y[0,i,0] = -(y1 - r*np.sin(theta1 + (theta2 - theta1)/(nu-1)*i)) # --- why??
        point_z[0,i,0] = 0
        xlow[i,0] = point_x[0,i,0]
        ylow[i,0] = point_y[0,i,0]
        
        
    # Top defn'
    for i in range(0, nu):
        point_x[nv-1,i,0] = x1 - r*np.cos(theta1 + (theta2 - theta1)/(nu - 1) * i )
        point_y[nv-1,i,0] = 2.0 + y1 - r*np.sin(theta1 + (theta2 - theta1)/(nu-1)* i)
        point_z[nv-1,i,0] = 0
        xhigh[i,0] = point_x[nv-1,i,0]
        yhigh[i,0] = point_y[nv-1,i,0]
        
    # For loop to implement Coon's patch! (thanks Caleb!)
    for k in range(0, nw):
        for j in range(0, nv):
            for i in range(0, nu):

                # What are these parameters?
                si = (1.0 * i ) / (nu - 1.0)
                sj = (1.0 * j ) / (nv - 1.0)

                # Cycling through x, y
                first_parta = (1.0 - si) * point_x[j,0,0] + sj * point_x[nv - 1, i, 0] + si * point_x[j, nu - 1,0] + (1 - sj) * point_x[0, i,0]
                second_parta = -(1.0 - sj) * (1.0 - si) * point_x[0,0,0] - sj * (1 - si) * point_x[nv - 1, 0,0] - sj * si * point_x[nv - 1, nu - 1,0] - (1-sj) * si * point_x[0,nu-1,0]
                point_x[j,i,0] = first_parta + second_parta

                first_partb = (1.0 - si) * point_y[j,0,0] + sj * point_y[nv-1,i,0] + si * point_y[j,nu-1,0] + (1-sj) * point_y[0,i,0]
                second_partb = -(1.0 - si) * (1 - si) * point_y[0,0,0] - sj * (1 - si) * point_y[nv-1,0,0] - sj * si * point_y[nv-1,nu-1,0] - (1 - sj) * si * point_y[0,nu-1,0]
                point_y[j,i,0] = first_partb + second_partb

                first_partc = (1.0 - si) * point_z[j,0,0] + sj * point_z[nv-1,i,0] + si * point_z[j,nu - 1,0] + (1 - sj) * point_z[0,i,0]
                second_partc = -(1.0 - sj) * (1 - si) * point_z[0,0,0] - sj * (1 - si) * point_z[nv-1,0,0] - sj * si * point_z[nv -1, nu - 1, 0] - (1 - sj) * si * point_z[0,nu - 1,0]
                point_z[j,i,k] = first_partc + second_partc


    # Converting from points to VTK readable format!
    #gridToVTK("./structured", point_x, point_y, point_z)

    # Compute the areas!
    areas = compute_areas(point_x, point_y)

    return point_x, point_y, point_z, xlow, ylow, xhigh, yhigh, areas