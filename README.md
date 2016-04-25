------------------------------------------------------------------
  Tungsten (ver. 1.0.0)
------------------------------------------------------------------

2D FINITE VOLUME SOLVER IN PYTHON

Currently, the code has the following characteristics:
   - Lax-Friedrichs method
   - Time derivatives via central differencing
   - Stable solution via artificial viscosity
   - 2D structured grid generation

Code uses evtk.hl libraries for converting structured data to a paraview file for post-processing.

Pranay Seshadri <br>
University of Cambridge <br>
ps583 <at> cam.ac.uk
<br>
To do:<br>
  - Convergence
  - Debugging
