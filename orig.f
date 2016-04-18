C*****************************************************************************C
C                                                                             C
C      THIS PROGRAM SOLVES THE EULER EQUATIONS FOR TWO_DIMENSIONAL DUCT FLOW  C
C      USING AN 'H' MESH WITH CELL CORNER STORAGE OF THE VARIABLES.           C
C                                                                             C
C*****************************************************************************C
C
C  First include the common blocks needed to transfer all variables.
C  and arrays between the subroutines.
C
      include 'common_block_euler'
      REAL FRKUT
C
C     You can plot the values in any array at any time by calling
C     subroutine "OUTPUT" with the required array as an argument. e.g.
C     CALL OUTPUT(AREA) to plot the array AREA . Then quit "EULER"
C     and run the program "EULPLT" .  EULPLT enables you to plot out the
C     the current values of all flow variables and Item 9 on the menu
C     of "EULPLT" plots the array you set as the argument of "OUTPUT".
C
C
C  Open the file "euler.log" which is used to plot the convergence
C  of the calculation via the separate program "PLTCONV".

      OPEN(UNIT=3,FILE='euler.log')
C
C  Call subroutine "read_data" to read in the data on the duct
C  geometry and flow conditions.
C
      call read_data
C      STOP
C
C    Call subroutine "inguess" this gives some arbitrary values to the
C    flow variables so that you can plot the grid. These values are not used
C    at all in the subsequent calculations.
C
      call inguess
C
C    Call subroutine "generate_grid" to set up the grid coordinates,
C    element areas and projected lengths of the sides of the elements.
C
      call generate_grid

c      call output(area)
c       STOP
C
C     You can call subroutine "OUTPUT" here to plot out the grid
C     you have generated.
C
c      call output(area)
c      STOP
C
C     Call subroutine "check_grid" to check that the areas and projected
C     lengths are correct.
C
      call check_grid
C
C     Subroutine  "crude_guess" is what its name says. It enables you to
C     start a calculation and obtain a solution but it will take longer than
C     necessary. When your program is working you should replace it
C     with "flow_guess" to obtain a better guess and a faster solution.
C
C
c      call crude_guess
C
      call new_guess
C
C      You can call "output" here to plot out your initial guess of
C      the flow field.
C
c      call output(ro)
c      stop
C
C     Call subroutine "set_timestep" to set the length of the timestep.
C     Initially this is a constant time step based on a conservative guess
C     of the Mach number.
C
      call set_timestep
C     STOP
C
C
C************************************************************************
C      START THE TIME STEPPING DO LOOP FOR "NSTEPS" LOOPS.
C************************************************************************
C
c     call set_others
c      call apply_bconds
      DO 1000 NSTEP = 1,NSTEPS
      DO I = 1, NI
       DO J = 1, NJ  
	 RO_START(I,J) = RO(I,J)
	 ROVX_START(I,J) = ROVX(I,J)
         ROVY_START(I,J) = ROVY(I,J)
         ROE_START(I,J) = ROE(I,J)
       END DO 
      END DO
      DO 999 NRKUT=1,4
      FRKUT = 1.0/(1 + 4 - NRKUT)
C      WRITE(*,*) FRKUT, NRKUT,RO(1,1)
C
      call set_others
      call apply_bconds
C
      call set_fluxes
C
      call sum_fluxes(frkut,fluxi_mass,fluxj_mass,ro,ro_start,delro)
      call sum_fluxes(frkut,fluxi_enth,fluxj_enth,roe,roe_start,delroe)
      call sum_fluxes(frkut,fluxi_xmom,fluxj_xmom,rovx,rovx_start,
     & delrovx)
      call sum_fluxes(frkut,fluxi_ymom,fluxj_ymom,rovy,rovy_start,
     & delrovy)
  999 CONTINUE
      call smooth(ro, corr_ro)
      call smooth(rovx, corr_rovx)
      call smooth(rovy, corr_rovy)
      call smooth(roe, corr_roe)
C
C

     
  
C
C     CHECK CONVERGENCE AND WRITE OUT SUMMARY EVERY 5 STEPS
C
      IF(MOD(NSTEP,5).EQ.0) THEN
C
	 call check_conv
	 call set_timestep
C
      ENDIF
C
C     GO TO 2000 IF CONVERGED TO THE INPUT TOLERANCE  "conlim"
C
      IF(EMAX.LT.CONLIM.AND.EAVG.LT.(0.5*CONLIM)) GO TO 2000
C
C     END OF THE MAIN TIME STEPPING DO LOOP
C
 1000 CONTINUE
C
C******************************************************************************
C
 2000 CONTINUE
C
      IF(EMAX.LT.CONLIM.AND.EAVG.LT.(0.5*CONLIM)) THEN
	  WRITE(6,*) ' CALCULATION CONVERGED IN ',NSTEP,' ITERATIONS'
	  WRITE(6,*) ' TO A CONVERGENCE LIMIT OF ', CONLIM
      ENDIF
C
C   CALCULATION FINISHED. CALL "output" TO WRITE THE PLOTTING FILE.
C
c      call output(diffrovx)
      call output(diffrovx)
C
      stop
      end
C
C******************************************************************************
C******************************************************************************
C
      subroutine read_data
C
      include 'common_block_euler'
C
C  Assign Fortran unit 1 to the file 'geom'
C  Assign Fortran unit 2 to file 'flow'
      OPEN(UNIT=1, FILE='test6_geom')
      OPEN(UNIT=2, FILE='test6_flow')
C
C
C
C      INSET YOUR CODE HERE
C
C  Read in the title and ni and nj from unit 1.
C  Before you do this check the format of the data file "geom".
C
c
      READ(1,*)  TITLE
      READ(1,*) NI,NJ
      
C     INSERT YOUR CODE HERE TO READ IN  ni and nj.
c
      IF(NI.LE.I_MAX.AND.NJ.LE.J_MAX)THEN
          DO I=1,NI
             READ(1,*) XLOW(I),YLOW(I),XHIGH(I),YHIGH(I)
C            WRITE(*,*) I , XLOW(I), YLOW(I), XHIGH(I), YHIGH(I)
          END DO
      ENDIF
C
C   Check that ni and nj are less than the array dimensions i_max and j_max.
c   Then read in xlow, ylow, xhigh, yhigh for each i  between 1 and ni
c
C      INSERT YOUR CODE HERE
C
C     Now read in the flow data from unit 2.
C     You should read in the following variables sequentially
      READ(2,*)  RGAS,GAMMA
      READ(2,*) PSTAGIN,TSTAGIN,ALPHA1,PDOWN
      READ(2,*) CFL,SMOOTH_FAC_IN
      READ(2,*) NSTEPS,CONLIM_IN
C
C
C     Set some other variables that will be used throughout the
C     calculation. Change ALPHA1 to radians. Scale the smoothing factor
C     and convergence limits by CFL.
C
      EMAX       = 1000000.
      EAVG       = EMAX
      CP         = RGAS*GAMMA/(GAMMA-1.)
      CV         = CP/GAMMA
      FGA        = (GAMMA - 1.)/GAMMA
      SMOOTH_FAC = SMOOTH_FAC_IN*CFL
      CONLIM     = CONLIM_IN*CFL
      ALPHA1     = ALPHA1*3.14159/180.
C
      CLOSE(1)
      CLOSE(2)
C
C      WRITE(*,*) YHIGH
      RETURN
      END
C
C******************************************************************************
C
      subroutine inguess
C
C  You do not need to touch this subroutine.
C
C
      include 'common_block_euler'
C
C  Set some completely arbitary values to enough of the flow
C  variables to enable plotting of the grid. These are not used at all
C  in the later calculations.
C
      DO I=1,NI
      DO J=1,NJ
	   RO(I,J)    = 1.2
	   ROVX(I,J)  = 100.*FLOAT(I)/NI
	   ROVY(I,J)  = 0.0
	   P(I,J)     = 100000.*(0.9 + 0.1*FLOAT(I)/NI)
	   HSTAG(I,J) = 300000.
	   ROE(I,J)   = P(I,J)/(GAMMA-1.)
      END DO
      END DO
C
C     CALCULATE THE REFERENCE VALUES WHICH ARE USED TO CHECK CONVERGENCE
C
      NCELLS    = NI*NJ
      JMID      = (1 + NJ)/2
      ROIN      = PSTAGIN/RGAS/TSTAGIN
      REF_RO    = (PSTAGIN-PDOWN)/RGAS/TSTAGIN
      REF_T     = TSTAGIN*(PDOWN/PSTAGIN)**FGA
      REF_V     = SQRT(2*CP*(TSTAGIN-REF_T))
      REF_ROVX   = ROIN*REF_V
      REF_ROVY   = REF_ROVX
      REF_ROE   = ROIN*CV*(TSTAGIN-REF_T)

      RETURN
      END
C
C******************************************************************************
C
      subroutine generate_grid
C
      include 'common_block_euler'
C
C   Calculate x and y values  x(i,j),y(i,j) of the grid nodes.
C
C   For each value of "i" the i-grid line joins (xlow(i),ylow(i)) to
C   (xhigh(i),yhigh(i)). For each value of "i" grid points (nodes) should be
C   linearly interpolated between these values for j between 1 and nj.
C   i.e.  x(i,1) should be xlow(i), x(i,nj) should be xhigh(i), etc.
c   
C   declare some intermediate variables
C
      REAL XO,YO,DEL_Y,DEL_X,DMI(I_MAX,J_MAX),DMJ(I_MAX,J_MAX)
c      WRITE(*,*) NI,NJ
       DO I = 1,NI
        DEL_Y = YHIGH(I) - YLOW(I)
        DEL_X = XHIGH(I) - XLOW(I)
C
        XO = XLOW(I)
        YO = YLOW(I) 
C      
        DO J = 1,NJ
            X(I,J) = XO + DEL_X * (J-1)/(NJ-1) 
            Y(I,J) = YO + DEL_Y * (J-1)/(NJ-1)
        END DO
      END DO


C     INSERT YOUR CODE HERE
C
c
      DO I = 1,NI-1
        DO J = 1,NJ-1
        A1 =  X(I,J+1)-X(I+1,J)
        B1 =  Y(I,J+1)-Y(I+1,J)
        A2 =  X(I+1,J+1)-X(I,J)
        B2 =  Y(I+1,J+1)-Y(I,J)
        C1 = ABS(A1 * B2)
        C2 = A2 * B1
        AREA(I,J) =  0.5*(C1 + C2)
c           WRITE(*,*) C1,C2,  AREA(I,J)
        END DO
      END DO

c    Calculate the areas of the cells AREA(I,J)
C    (N.B. There are (ni-1) x (nj-1) cells.
c  The area of a quadrilateral (regular or irregular) can be shown to be
c  half of the cross product of the vectors forming the diagonals.
c  See Hirsch volume 1, section 6.2.1. (or lecture).
c  Make sure that the area comes out positive !
C
C    CODE IS MODIFIED SO DMIN IS NOW NOT A SINGLE VALUE BUT A MATRIX (I,J)
      DMIN_I = 100
      DMIN_J = 100
       DO I = 1,NI
          DO J = 1,NJ-1
           DLIX(I,J) = Y(I,J+1) - Y(I,J)
           DLIY(I,J) = -X(I,J+1) + X(I,J)
           DMI(I,J)= (DLIX(I,J)**2+DLIY(I,J)**2)**0.5
          END DO
       END DO 

       DO I = 1,NI-1
         DO J = 1,NJ
          DLJX(I,J) = -Y(I+1,J) + Y(I,J)
          DLJY(I,J) =  X(I+1,J) - X(I,J)
 	   DMJ(I,J)= (DLJX(I,J)**2+DLJY(I,J)**2)**0.5
         END DO
       END DO
C      
       DO I=1,NI-1
        DO J=1,NJ-1
        DMIN(I,J) = MIN(DMI(I,J),DMJ(I,J))
C	WRITE(*,*) DMI(I,J) , DMJ(I,J)
       END DO
       END DO
c      STOP
C
       DO I=1,NI-1
          DMIN(I,NJ) = DMJ(I,NJ)
       END DO
C      
       DO J=1,NJ-1
          DMIN(NI,J) = DMI(NI,J)
       END DO
C      
C      
        DMIN(NI,NJ) = MIN( DMJ(NI-1,NJ-1),DMI(NI-1,NJ-1))
c       WRITE(*,*) DMIN
c   Calculate the x and y components of the length vector of the i-faces
c   (i.e. those corresponding to i = constant).
c   The length vector of a face is a vector normal to the face with
C   magnitude equal to the length of the face.
C   It is positive in the direction of an inward normal to the cell i,j .
c   Call these lengths dlix(i,j) and dliy(i,j)
c
C    INSERT YOUR CODE HERE
C
c   Now calculate the x and y components of the length vector of the j-faces.
c   (i.e. those corresponding to j = constant)
c   Call these lengths dljx(i,j) and dljy(i,j)
C
C    INSERT YOUR CODE HERE
C
C    Now find "DMIN" the minimum length scale of any element. This is
C    defined as the length of the shortest side of the element.
C    Call this minimum "DMIN". It is used to set the time step.
C
C     INSERT YOUR CODE HERE
C

C      write(*,*) XLOW
c
C
      return
      end
C
C******************************************************************************
C******************************************************************************
C
      subroutine check_grid
C
      include 'common_block_euler'
C	
      REAL SMVAL,SUMMATION_1,SUMMATION_2,TVAL
C
C   Check your grid and areas for both the "bump" and the "bend"
C   test data.
C
C  First check that all areas are positive (by program or writing out)
C
C     INSERT YOUR CODE HERE
      SMVAL = -0.0001
      TVAL = 0.0001
      DO I=1,NI-1
	 DO J=1,NJ-1
	   IF(AREA(I,J).LT.SMVAL)THEN
		 WRITE(*,*) 'ERROR AREA IS NEGATIVE'
 	  ENDIF
 	END DO
      END DO

C   
C
C  Next check that the sum of the length vectors of the 4 faces
C  of every element is very nearly zero in each coordinate direction.
C  It is absolutely essential that this is correct !
C
      DO I = 1,NI-1
	 DO J = 1,NJ-1
           SUMMATION_1 = DLIX(I,J)-DLIX(I+1,J)+DLJX(I,J)- DLJX(I,J+1)
           SUMMATION_2 = DLIY(I,J)-DLIY(I+1,J)+DLJY(I,J)- DLJY(I,J+1)
	   IF(SUMMATION_1.GT.TVAL.OR.SUMMATION_2.GT.TVAL)THEN
		WRITE(*,*) 'LENGTH OF 4 VECTORS IS NOT ZERO!'
           ENDIF
	 END DO
       END DO	
c      STOP
C  If not go back and check your subroutine "generate_grid".
C
C  Be careful with a test of the form
C             if( a.eq.0.0 ) then .....
C  This will probably never be true.  Computers work to a finite number of
C  significant figures and "a" will probably be +0.0000001 or -0.0000001.
C  Test for something like
C             if( abs(a).le.small_number ) then ...
C
C
C      INSERT YOUR CODE HERE
C
C  Any other tests that you can think of. For example you could plot
C  contours of  "area(i,j)" by using  --  call output(area) .
C
C
      return
      end                                                                   
C                                             
C******************************************************************************
C******************************************************************************
C
      subroutine flow_guess
C
      include 'common_block_euler'
C
      DIMENSION AFLOW(I_MAX), V_GUESS(I_MAX), RO_GUESS(I_MAX)
      REAL    MFLOW,MLIM
C
C   In this subroutine we make an initial guess of the primary variables
C   i.e.  RO, ROVX, ROVY and ROE.
C   The guess does not need to be very accurate but the better
C   it is the faster your program will converge.
C   You should assign values to RO(I,J), ROVX(I,J), ROVY(I,J) and ROE(I,J)
C   at every grid point in this subroutine.
C
C     Work out the length of each "i" line between grid points
C     "i,1" and "i,nj" and call it  "AFLOW(I)" .
C
C      INSERT YOUR CODE HERE

	DO I=1,NI
	AFLOW(I) = SQRT((XHIGH(i)-XLOW(i))**2+(YHIGH(i)-YLOW(i))**2)
	enddo
C
C   Make an initial guess of the density and velocity at the exit by
C   assuming isentropic flow between the inlet stagnation pressure PSTAGIN
C   and temperature TSTAGIN and the exit static pressure PDOWN.
C   Use these together with "AFLOW(NI)" to estimate the mass flow rate.
C   Call this "MFLOW".
C
C      INSERT YOUR CODE HERE

	ROSTAGIN = PSTAGIN/(RGAS*TSTAGIN)
	ROEXIT = ROSTAGIN*(PDOWN/PSTAGIN)**(1./GAMMA)
	TEXIT = PDOWN/(RGAS*ROEXIT)
	VEXIT = SQRT(2*CP*(TSTAGIN-TEXIT))

	MFLOW = ROEXIT*AFLOW(NI)*VEXIT

	
C
C
C     SET A LIMIT TO THE MAXIMUM ALLOWABLE MACH NUMBER IN THE INITIAL
C     GUESS. CALL THIS "MACHLIM". CALCULATE THE CORRESPONDING TEMPERATURE.
C
      MACHLIM = 1.0
      TLIM = TSTAGIN/(1.0 + 0.5*(GAMMA-1.0)*MACHLIM*MACHLIM)
C
C    Now estimate the velocity and density at every "i" line.
C    Call the velocity V_GUESS(I) and the density RO_GUESS(I).
C
C    First assume that the density is constant and equal to the exit
C    density calculated above and that the flow is perpendicular to the
C    "i" = constant lines and hence occupies area AFLOW(I).
C    Hence use continuity to estimate the flow velocity V_GUESS(I).
C    Use this velocity to calculate the static temperature assuming
C    that the stagnation temperature is constant.
C    Check that this temperature is not less than TLIM and set = TLIM
C    if it is.
C    Next use this temperature and isentropic flow to obtain a better
C    estimate of the density, RO_GUESS(I).
C    Use this density and continuity to obtain a better estimate of
C    the velocity, set = V_GUESS(I).
C
C      INSERT YOUR CODE HERE

	DO I=1,NI-1
		
	RO_GUESS(i) = ROEXIT
	V_GUESS(i) = MFLOW/(RO_GUESS(i)*AFLOW(i))
	TSTATIC = TSTAGIN - 0.5*((V_GUESS(i)**2)/CP)

	IF (TSTATIC.LT.TLIM) THEN
	TSTATIC = TLIM
	ENDIF

	RO_GUESS(i) = ROSTAGIN*(TSTATIC/TSTAGIN)**(1./(GAMMA-1))
	V_GUESS(i) = MFLOW/(RO_GUESS(i)*AFLOW(i))


	END DO	

C
C   Direct the velocity found above along the J= constant grid lines to find
C   the velocity VX(I,J) in the  x  direction and VY(I,J) in the  y  direction.
C   Use these and RO_GUESS(I) to set ROVX(I,J), ROVY(I,J) and ROE(I,J).
C   Also set RO(I,J).
C   Note that ROE(I,J) includes the kinetic energy component of the
C   internal energy.
C
C     INSERT YOUR CODE HERE

	
	DO J = 1,NJ
	DO I = 1,NI-1
	

	DX = X(i+1,j) - X(i,j)
	DY = Y(i+1,j) - Y(i,j)
	DS = SQRT(DX*DX+DY*DY)

	VX(i,j) = V_GUESS(i)*DX/DS
	VY(i,j) = V_GUESS(i)*DY/DS
	ROVX(i,j) = VX(i,j)*RO_GUESS(i)
	ROVY(i,j) = VY(i,j)*RO_GUESS(i)

	NUME = 0.5*(VX(i,j)*VX(i,j)+VY(i,j)*VY(i,j))
	TSTATIC = TSTAGIN - NUME/Cp
	ROE(i,j) = RO_GUESS(i)*(CV*TSTATIC + NUME)
	RO(i,j) = RO_GUESS(i)	


	enddo

      	ROVX(NI,J) = ROVX(NI-1,J)
	ROVY(NI,J) = ROVY(NI-1,J)
      	RO(NI,J)   = RO(NI-1,J)
      	ROE(NI,J)  = ROE(NI-1,J)
	
	enddo

C
C  Store the "old" values of the variables for use in the first
C  convergence check in subroutine "check_conv"
C
      DO I=1,NI
      DO J=1,NJ
	      RO_OLD(I,J)   = RO(I,J)
	      ROVX_OLD(I,J) = ROVX(I,J)
	      ROVY_OLD(I,J) = ROVY(I,J)
	      ROE_OLD(I,J)  = ROE(I,J)
      END DO
      END DO
C
      RETURN
      END
C

C******************************************************************************
C******************************************************************************
C
      SUBROUTINE CRUDE_GUESS
C     You should not need to touch this subroutine
C
C     This subroutine makes a very simple and inaccurate guess of the
C     flow field. Enough to get the program working.
C
      INCLUDE 'common_block_euler'
C
      JMID   = NJ/2
      TDOWN  = TSTAGIN*(PDOWN/PSTAGIN)**FGA
      VDOWN  = SQRT(2*CP*(TSTAGIN - TDOWN))
      RODOWN = PDOWN/RGAS/TDOWN
C
      DO 20 J=1,NJ
      DO 10 I=1,NI-1
      DX  = X(I+1,JMID) - X(I,JMID)
      DY  = Y(I+1,JMID) - Y(I,JMID)
      DS  = SQRT(DX*DX + DY*DY)
      XVEL      = VDOWN*DX/DS
      YVEL      = VDOWN*DY/DS
      ROVX(I,J) = RODOWN*XVEL
      ROVY(I,J) = RODOWN*YVEL
      RO(I,J)   = RODOWN
      ROE(I,J)  = RODOWN*(CV*TDOWN + 0.5*VDOWN*VDOWN)
   10 CONTINUE
      ROVX(NI,J) = ROVX(NI-1,J)
      ROVY(NI,J) = ROVY(NI-1,J)
      RO(NI,J)   = RO(NI-1,J)
      ROE(NI,J)  = ROE(NI-1,J)
   20 CONTINUE
      RETURN
      END
C
C**********************************************************************
C**********************************************************************
C
      subroutine set_others
C
      include 'common_block_euler'
C
C  This routine calculates secondary flow variables from the primary ones
C  at every grid point.
C
C   The primary variables are RO,ROVX,ROVY and ROE
C
C  The secondary variables are the velocity components VX(I,J) and VY(I,J),
C  the static pressure P(I,J) and the stagnation enthalpy HSTAG(I,J).
C  Note:  "HSTAG"  NOT  "HO".
C

       DO I=1,NI          
        DO J=1,NJ
         VX(I,J) = ROVX(I,J)/RO(I,J)
         VY(I,J) = ROVY(I,J)/RO(I,J)
c        CLASS C MODIFICATION -- SINCE ROE IS NOW NO LONGER BEING CALCULATED SO 
C        CANNOT BE USED TO CALCULATE ANY OTHER VARIABLES
c	 E = ROE(I,J)/RO(I,J)  
c	 T = (E-0.5*VX(I,J)**2-0.5*VY(I,J)**2)/CV
c	 T = TO * (RO(I,J) / (PO/(RGAS * TO) ) )**(GAMMA - 1_)
C         WRITE(*,*) T
C         HO = CP * TSTAGIN
C
         HSTAG(I,J) = CP*TSTAGIN
	 KE = 0.5* (VX(I,J)*VX(I,J) + VY(I,J)*VY(I,J))
        T = (HSTAG(I,J) - KE)/(CP)
         P(I,J) = RO(I,J)*RGAS*T
        END DO
       END DO

       RETURN
       END
C******************************************************************************
C******************************************************************************

      subroutine set_timestep
C
      include 'common_block_euler'
      REAL VELO(I_MAX,J_MAX),VEL_AVG,A_AVG
C
C  This subroutine sets the length of the time step based on the
C  stagnation speed of sound "ASTAG" and the minimum length scale
C  of any element, "DMIN". The timestep must be called "DELTAT"
C
C  An assumption that the maximum flow speed will be equal to "ASTAG"
C  is also made. This will be pessimistic for subsonic flows
C  but may be optimistic for supersonic flows. In the latter case the
C  length of the time step as determined by "CFL" may need to be reduced.
C
C  The CFL number was input as data in data set "FLOW"
C
c     CLASS C MODIFICATION ----------------------------
       DO I = 1,NI
       DO J = 1,NJ
          VELO(I,J) =SQRT(VX(I,J)*VX(I,J) + VY(I,J)*VY(I,J))
       END DO
       END DO
C
      DO I=1,NI-1         
      DO J=1,NJ-1
C        LOCAL SPEED OF SOUND IS THE AVERAGE OF THE FOUR POINTS: 
         A_1 = SQRT(GAMMA * P(I,J)/RO(I,J))
         A_2 = SQRT(GAMMA * P(I+1,J)/RO(I+1,J) )
         A_3 = SQRT(GAMMA * P(I,J+1)/RO(I,J+1) )
         A_4 = SQRT(GAMMA * P(I+1,J+1)/RO(I+1,J+1))
         A_AVG = 0.25*(A_1 + A_2 + A_3 + A_4)
C
C        LOCAL VELOCITY IS THE AVERAGE OF THE FOUR POINTS:
         VEL_1 = VELO(I,J)
         VEL_2 = VELO(I+1,J)
         VEL_3 = VELO(I,J+1)
         VEL_4 = VELO(I+1,J+1)
         VEL_AVG = 0.25*(VEL_1 + VEL_2 + VEL_3 + VEL_4) 
c         WRITE(*,*) DMIN(I,J), A_AVG, VEL_AVG
         STEP(I,J) =  CFL * DMIN(I,J)/(A_AVG + VEL_AVG)
C	  WRITE(*,*) DMIN(I,J), A_AVG, VEL_AVG
      END DO
      END DO
C      STOP
C
C     NOW WE HAVE TO MANUALLY INPUT VALUES FOR STEP(I,J) AT THE LAST ROW 
C     AND COLUMN
      DO I=1,NI-1
         A_AVG = SQRT(GAMMA * P(I,NJ)/RO(I,NJ))
         VEL_AVG = VELO(I,NJ)
         STEP(I,NJ) = CFL * DMIN(I,NJ)/(A_AVG + VEL_AVG)
      END DO
C      
      DO J=1,NJ-1
         A_AVG = SQRT(GAMMA * P(NI,J)/RO(NI,J))
         VEL_AVG = VELO(NI,J)
         STEP(NI,J) = CFL * DMIN(NI,J)/(A_AVG + VEL_AVG)
      END DO
C      
      A_AVG = SQRT(GAMMA * P(NI,NJ)/RO(NI,NJ))
      STEP(NI,NJ) = CFL * DMIN(NI,J)/(VELO(NI,NJ)+A_AVG)
C
	
	
C      DELTAT = 1
      RETURN
      END
C      
      


C*****************************************************************************
C
      subroutine set_fluxes
C
      include 'common_block_euler'
C
C  This subroutine calculates the fluxes of mass momentum and energy
C  across the faces of every cell.
C
C  The "i" faces with i = 1 or i = ni are the upstream and downstream
C  boundaries to the flow domain, while "j" faces with j = 1 or j = nj
C  are solid boundaries. All fluxes are calculated assuming a linear
C  variation in the flow properties between the cell corner nodes.
C
C
C  First calculate the mass flux across each "i" face of the elements.
C  Also calculate the total mass flow rate "flow(i)" across each "i" line.
C  This will be used to check global continuity.
C
      do i=1,ni
      flow(i) = 0.0
      do j=1,nj-1
	fluxi_mass(i,j) =  0.5*(  (rovx(i,j)+rovx(i,j+1))*dlix(i,j) +
     &                            (rovy(i,j)+rovy(i,j+1))*dliy(i,j)  )

	flow(i) = flow(i) + fluxi_mass(i,j)

      end do
      end do

c
c     Now the mass flux across each "j" face.
c
      do i=1,ni-1
      do j=2,nj
C
C           INSERT YOUR CODE HERE TO CALCULATE "fluxj_mass(i,j)"
C
        fluxj_mass(i,j) =  0.5*(  (rovx(i,j)+rovx(i+1,j))*dljx(i,j) +
     &                            (rovy(i,j)+rovy(i+1,j))*dljy(i,j)  )
      end do
      end do
C
C   Set the mass fluxes through the j=1 and j=nj faces to zero as
C   these are solid surfaces. It is not necessary to resolve the
C   velocity parallel to the surfaces.
C
      do i=1,ni-1
         fluxj_mass(i,1) = 0.0
         fluxj_mass(i,nj)= 0.0
      end do
c
c       Calculate the fluxes of X-momentum---------------------
c
      do i=1,ni
         do j=1,nj-1
	fluxi_xmom(i,j) = 0.5*(   fluxi_mass(i,j)*(vx(i,j)+vx(i,j+1))  +
     &                            (p(i,j)+p(i,j+1))*dlix(i,j)   )
          end do
      end do
c
     
C
C      INSERT YOUR CODE HERE TO SET "fluxj_xmom(i,j)"
C
      do i=1,ni-1
        do j=1,nj
        fluxj_xmom(i,j) = 0.5*(   fluxj_mass(i,j)*(vx(i,j)+vx(i+1,j))  +
     &                            (p(i,j)+p(i+1,j))*dljx(i,j)   )
         end do
      end do
C
c       Calculate the fluxes of Y-momentum---------------------
C
C      INSERT YOUR CODE HERE TO SET "fluxi_ymom(i,j)"
       do i=1,ni
         do j=1,nj-1
	fluxi_ymom(i,j) = 0.5*(   fluxi_mass(i,j)*(vy(i,j)+vy(i,j+1))  +
     &                            (p(i,j)+p(i,j+1))*dliy(i,j)   )
          end do
      end do
C
C      INSERT YOUR CODE HERE TO SET "fluxj_ymom(i,j)"
       do i=1,ni-1
        do j=1,nj
        fluxj_ymom(i,j) = 0.5*(   fluxj_mass(i,j)*(vy(i,j)+vy(i+1,j))  +
     &                            (p(i,j)+p(i+1,j))*dljy(i,j)   )
         end do
      end do
c
c       Calculate the fluxes of Enthalpy-----------------------
c
C
C      INSERT YOUR CODE HERE TO SET "fluxi_enth(i,j)"
      do i=1,ni
         do j=1,nj-1
	fluxi_enth(i,j) = 0.5*fluxi_mass(i,j)*(HSTAG(i,j)+HSTAG(i,j+1))
          end do
      end do

C
C      INSERT YOUR CODE HERE TO SET "fluxj_enth(i,j)"
      do i=1,ni-1
         do j=1,nj
	fluxj_enth(i,j) = 0.5*fluxj_mass(i,j)*(HSTAG(i,j)+HSTAG(i+1,j))
          end do
      end do
C
C   Note that we could have set the flux of enthalpy to zero on
C   j=1 and j=nj. This would save a bit of CPU time but the fluxes
C   will be zero anyhow since the mass fluxes were set to zero.
C
      return
      end
C
C**************************************************************************
C**************************************************************************
C
      subroutine sum_fluxes (frkut,iflux,jflux,prop,prop_start,delprop)
C
C   This subroutine sums the fluxes for each element, calculates
C   the changes in the variable "PROP" and distributes them to the
C   four corners of the element.
C
      INCLUDE 'common_block_euler'
C
      REAL       IFLUX(I_MAX,J_MAX),JFLUX(I_MAX,J_MAX),
     &           PROP(I_MAX,J_MAX),DELPROP(I_MAX,J_MAX),
     &           STORE(I_MAX,J_MAX),PROP_START(I_MAX,J_MAX)

C
C   Find the change in the variable "PROP" in each cell over the
C   time step "DELTAT" and save it in "STORE(I,J)"
C     
      DO I=1,NI-1
      DO J=1,NJ-1
C     
      TOTFLUX=(IFLUX(I,J) - IFLUX(I+1,J) + JFLUX(I,J) - JFLUX(I,J+1))
      STORE(I,J) = FRKUT*STEP(I,J)* TOTFLUX/area(i,j)
      

C     INSERT YOUR CODE HERE to calculate the change in the variable
C     "PROP" over the time step "deltat" and sets the result to "STORE(I,J)"
C
      END DO
      END DO
C
C  Now distribute the changes equally to the four interior corners of each
C  cell. Each interior grid points receive one quarter of the change
C  from each of the four cells adjacent to it.
C
      DO I=2,NI-1
      DO J=2,NJ-1
C
C     INSERT YOUR CODE HERE to calculate "ADD" the change to be
C     added to each interior node.
      ADD = (STORE(I,J)+STORE(I-1,J-1)+STORE(I,J-1)+STORE(I-1,J))*0.25
      PROP(I,J)     = PROP_START(I,J)   + FRKUT*ADD
      END DO
      END DO
C
C  Now add the changes to the upper and lower boundaries.
C  These receive half the change from each of the two cells
C  adjacent to them.
C
      DO I=2,NI-1
C
C   INSERT YOUR CODE HERE to calculate "ADD" for the nodes with J=1.
C     
      ADD = (STORE(I,1)+STORE(I-1,1))*0.5
      PROP(I,1)     = PROP_START(I,1)    +   FRKUT*ADD
C
C   INSERT YOUR CODE HERE to calculate "ADD" for the nodes with J=NJ.
C
      ADD = (STORE(I,NJ-1)+STORE(I-1,NJ-1))*0.5
      PROP(I,NJ)    = PROP_START(I,NJ)   +   FRKUT*ADD
C
      END DO
C
C  Now add on extra changes to the inlet & outlet boundary points
C  These receive half the change from each of the two cells
C  adjacent to them.
C
      DO J=2,NJ-1
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the nodes with I=NI.
C
      ADD = (STORE(NI-1,J)+STORE(NI-1,J-1))*0.5
      PROP(NI,J)     = PROP_START(NI,J)   + FRKUT*ADD
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the nodes with I=1.
C
      ADD = (STORE(1,J)+STORE(1,J-1))*0.5
      PROP(1,J)      = PROP_START(1,J)    + FRKUT*ADD
C
      END DO
C
C      Now add the changes on to the four corner points.
C      These receive the full change from the single cell of which
C      they form one corner.
C
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the node with I=1,J=1.
C
      ADD = STORE(1,1)
      PROP(1,1)       = PROP_START(1,1)  +  FRKUT*ADD
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the node with I=1,J=NJ.
C
      ADD = STORE(1,NJ-1)
      PROP(1,NJ)      = PROP_START(1,NJ)  +  FRKUT*ADD
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the node with I=NI,J=1.
C
      ADD = STORE(NI-1,1)
      PROP(NI,1)      = PROP_START(NI,1) +  FRKUT*ADD
C
C     INSERT YOUR CODE HERE to calculate "ADD" for the node with I=NI,J=NJ.
C
      ADD = STORE(NI-1,NJ-1)
      PROP(NI,NJ)     = PROP_START(NI,NJ) + FRKUT*ADD
C
C
C Now save the changes in the primary variables as "DELPROP".
C This will be used in the convergence check and also in future
C improvements to the scheme.
C
      DO  I=1,NI-1
      DO  J=1,NJ-1
      DELPROP(I,J) = STORE(I,J)
      END DO
      END DO
C
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE SMOOTH (PROP, CORRPROP)
C
C
C  This subroutine smooths the variable "PROP" (i.e. it adds the
C  artificial viscosity) by taking (1-SF) x the calculated value of
C  "PROP" + SF x (The average of the surrounding values of "PROP").
C   Where SF is the smoothing factor.
C
      INCLUDE 'common_block_euler'
C
      DIMENSION  PROP(I_MAX,J_MAX),STORE(I_MAX,J_MAX),
     &           CORRPROP(I_MAX,J_MAX)
      REAL FCORR
      FCORR = 0.95
C
C   To avoid using already smoothed values for smoothing other values
C   the smoothed values are initially stored in an array "STORE".
C
      SF   = SMOOTH_FAC
      SFM1 = 1.0 - SF
C
      DO I = 1,NI
	   IP1  = I+1
	   IF(I.EQ.NI) IP1 = NI
	   IM1  = I-1
	   IF(I.EQ.1) IM1 =1
      DO J = 2,NJ-1
	   AVG = 0.25*(PROP(IP1,J)+PROP(IM1,J)+PROP(I,J-1)+PROP(I,J+1))
           CORRNEW = FCORR*(PROP(I,J) - AVG)
           CORRPROP(I,J) = 0.99*CORRPROP(I,J) + 0.01*CORRNEW
C
C      ADD YOUR CODE HERE. i.e. STORE(I,J) =  ????
C
          STORE(I,J) = SFM1*PROP(I,J)+SF*(AVG+CORRPROP(I,J))
      END DO
C
C   On the surfaces j=1 and j=nj take the average as shown below.
C
      AVG1  = (PROP(IM1,1) + PROP(IP1,1) + 2.*PROP(I,2) - PROP(I,3))/3.
      CORRNEW = FCORR * (PROP(I, 1) - AVG1)
      CORRPROP(I,1) = 0.99*CORRPROP(I,1) + 0.01*CORRNEW

      AVGNJ = (PROP(IM1,NJ) + PROP(IP1,NJ) + 2.*PROP(I,NJ-1)
     &      -  PROP(I,NJ-2))/3.
      CORRNEW = FCORR * (PROP(I,NJ) - AVGNJ)
      CORRPROP(I,NJ) = 0.99*CORRPROP(I,NJ) + 0.01*CORRNEW
C
C    ADD YOUR CODE HERE TO SMOOTH THE SURFACE VALUES
C    i.e.  STORE(I,1) =  ??
C          STORE(I,NJ) = ??
C
      STORE(I,1) = SFM1*PROP(I,1)+SF*(AVG1+CORRPROP(I,1))
      STORE(I,NJ) = SFM1*PROP(I,NJ)+SF*(AVGNJ+CORRPROP(I,NJ))
      END DO
C
C   Reset the smoothed value to "PROP" before returning to the main program.
C
      DO I=1,NI
      DO J=1,NJ
	   PROP(I,J) = STORE(I,J)
      END DO
      END DO
C
C
      RETURN
      END
C
C
C*****************************************************************************
C******************************************************************************
C
      subroutine apply_bconds
C
      INCLUDE 'common_block_euler'
	
      REAL VEL
C
C  This subroutine applies the boundary conditions that P = PDOWN
C  at I = NI. At the inlet boundary the change in density is relaxed
C  to obtain ROINLET(J) which is then used to obtain the other properties
C  at inlet assuming isentropic flow from stagnation conditions "PSTAGIN"
C  and TSTAGIN" together with the specified inlet flow angle "ALPHA1".

C
C  Because the inlet condition may become unstable it is safer to
C  relax the changes in inlet density by a factor "RFIN"
C  Typically "RFIN" = 0.25 as set below. Reduce this if the inlet
C  becomes unstable.
C
C  It is also worth checking if "ROINLET" is greater than "ROSTAGIN"
C  and setting ROINLET to 0.9999*ROSTAGIN if it is.
C  This saves the program crashing during transients.
C
C
      RFIN     = 0.25
      RFIN1    = 1.0-RFIN
      ROSTAGIN = PSTAGIN/(RGAS*TSTAGIN)
C
      DO J=1,NJ
      IF(NSTEP.EQ.1) THEN
	   ROINLET(J) = RO(1,J)
      ELSE
	   ROINLET(J) = RFIN*RO(1,J) + RFIN1*ROINLET(J)
      ENDIF
C
      IF(ROINLET(J).GT.0.9999*ROSTAGIN) ROINLET(J)=0.9999*ROSTAGIN
C
C      INSERT YOUR CODE HERE to calculate P(1,J),ROVX(1,J),ROVY(1,J)
C      and ROE(1,J)  from ROINLET(J), PSTAGIN, TSTAGIN  and ALPHA1.
C      also set VX(1,J), VY(1,J) and HSTAG(1,J)
C     

      P(1,J) = PSTAGIN*( ROINLET(J)/ROSTAGIN )**GAMMA
	T = P(1,J)/(RGAS*ROINLET(J))
      VEL = ( 2.0*CP*(TSTAGIN - T))**0.5
      VX(1,J) = VEL * COS(ALPHA1)
      VY(1,J) = VEL * SIN(ALPHA1)
      ROVX(1,J) = ROINLET(J) * VX(1,J)
      ROVY(1,J) = ROINLET(J) * VY(1,J)
      ROE(1,J) = ROINLET(J)*(CV*T + 0.5*VEL**2)
      HSTAG(1,J) = CP * TSTAGIN 
      P(NI,J) = PDOWN
C   
      END DO
C
C
C   Set the pressure at the downstream boundary i=ni to the exit
C   static pressure "PDOWN" for all j values.
C
C      INSERT YOUR CODE HERE
C	
      
C
      RETURN
      END
C*****************************************************************************
C
      subroutine check_conv
C
C       You should not need to change this subroutine
C
      INCLUDE 'common_block_euler'
C
C  This subroutine checks the changes in all primary variables over
C  the last 5 steps.
C
C  It also writes a short output summary to unit 6 and a file for
C  plottimg the convergence history to unit 3.

      DELROMAX   = 0.0
      DELROVXMAX = 0.0
      DELROVYMAX = 0.0
      DELROEMAX  = 0.0
      DELROAVG   = 0.0
      DELROVXAVG = 0.0
      DELROVYAVG = 0.0
      DELROEAVG  = 0.0
      IMAX = 0
      JMAX = 0
C
C     "IMAX,JMAX" is the grid point where the change in ROVX is a max.
C
      DO I=1,NI
      DO J=1,NJ
C
	   DELTA = ABS(RO(I,J) - RO_OLD(I,J))
	   IF(DELTA.GT.DELROMAX) DELROMAX = DELTA
	   DELROAVG = DELROAVG + DELTA
C
	   DELTA = ABS(ROVX(I,J)-ROVX_OLD(I,J))
	   IF(DELTA.GT.DELROVXMAX)  THEN
	       DELROVXMAX = DELTA
	       IMAX = I
	       JMAX = J
	   ENDIF
	   DIFFROVX(I,J) = DELTA
	   DELROVXAVG = DELROVXAVG + DELTA
C
	   DELTA     = ABS(ROVY(I,J) - ROVY_OLD(I,J))
	   IF(DELTA.GT.DELROVYMAX) DELROVYMAX = DELTA
	   DELROVYAVG = DELROVYAVG + DELTA
C
	   DELTA = ABS(ROE(I,J) - ROE_OLD(I,J))
	   IF(DELTA.GT.DELROEMAX) DELROEMAX = DELTA
	   DELROEAVG = DELROEAVG + DELTA
C
      END DO
      END DO
C
C     CALCULATE THE AVERAGE CHANGES
C
      DELROAVG   =  DELROAVG/NCELLS/REF_RO
      DELROVXAVG = DELROVXAVG/NCELLS/REF_ROVX
      DELROVYAVG = DELROVYAVG/NCELLS/REF_ROVY
      DELROEAVG  = DELROEAVG/NCELLS/REF_ROE
      DELROVXMAX = DELROVXMAX/REF_ROVX
      DELROVYMAX = DELROVYMAX/REF_ROVY
C
      EMAX      = AMAX1(DELROVXMAX,DELROVYMAX)
      EAVG      = AMAX1(DELROVXAVG,DELROVYAVG)
C
C   Store the maximum change in ROVX as EMAX to be printed out.
C
C    Save the current values of the primary variables as PROP_OLD values
C    for use in the next convergence check.
C
      DO I=1,NI
      DO J=1,NJ
	   RO_OLD(I,J)   = RO(I,J)
	   ROVX_OLD(I,J) = ROVX(I,J)
	   ROVY_OLD(I,J) = ROVY(I,J)
	   ROE_OLD(I,J)  = ROE(I,J)
      END DO
      END DO
C
C  Write the average changes in the primary variables to unit 3
C  for use in the convergence plotting program "PLTCONV"
C
      WRITE(3,300) DELROAVG,DELROVXAVG,DELROVYAVG,DELROEAVG
  300 FORMAT(4E13.6)
C
C  Write a short output summary to unit 6 which will usually be the screen.
C
      FLOW_RATIO = FLOW(NI)/FLOW(1)
      WRITE(6,*) ' TIME STEP NUMBER ', NSTEP
      WRITE(6,5) emax,imax,jmax,eavg
    5 FORMAT(' EMAX= ',E10.3,' AT IMAX = ',I5,' JMAX= ',I5,' EAVG= ',
     & E10.3)
      WRITE(6,*) 'inlet flow= ',FLOW(1),' outlet to inlet flow ratio',
     & FLOW_RATIO
C
      RETURN
      END
C
C------------------------------------------------------------------
      SUBROUTINE NEW_GUESS
C
C     This subroutine makes a guess of the flow field assuming no
C     variation in the cross-flow (J) direction and assuming a
C     linear variation of flow properties from inlet to outlet.
C
      INCLUDE 'common_block_euler'
C
      JMID   = NJ/2
      TDOWN  = TSTAGIN*(PDOWN/PSTAGIN)**FGA
      VDOWN  = SQRT(2*CP*(TSTAGIN - TDOWN))
      RODOWN = PDOWN/RGAS/TDOWN
      PINLET = 55000.
      TINLET  = TSTAGIN*(PINLET/PSTAGIN)**FGA
      VINLET  = SQRT(2*CP*(TSTAGIN - TINLET))
      ROIN    = PINLET/RGAS/TINLET

C
      DO 20 J=1,NJ
      DO 10 I=1,NI-1
      DX  = X(I+1,JMID) - X(I,JMID)
      DY  = Y(I+1,JMID) - Y(I,JMID)
      DS  = SQRT(DX*DX + DY*DY)
      VLOCAL    = VINLET  + (VDOWN-VINLET)*FLOAT(I-1)/(NI-1)
      ROLOCAL   = ROIN    + (RODOWN-ROIN)*FLOAT(I-1)/(NI-1)
      TLOCAL    = TINLET  + (TDOWN-TINLET)*FLOAT(I-1)/(NI-1)

      XVEL      = VLOCAL*DX/DS
      YVEL      = VLOCAL*DY/DS
      ROVX(I,J) = ROLOCAL*XVEL
      ROVY(I,J) = ROLOCAL*YVEL
      RO(I,J)   = ROLOCAL
      ROE(I,J)  = ROLOCAL*(CV*TLOCAL + 0.5*VLOCAL*VLOCAL)
   10 CONTINUE
      ROVX(NI,J) = ROVX(NI-1,J)
      ROVY(NI,J) = ROVY(NI-1,J)
      RO(NI,J)  = RO(NI-1,J)
      ROE(NI,J) = ROE(NI-1,J)
   20 CONTINUE
      RETURN
      END

C******************************************************************************
C******************************************************************************
C
      SUBROUTINE OUTPUT(PLOTVAR)
C
C            You should not need to touch this subroutine.
C
C
C      This subroutine writes a file to unit 7 for use by the plotting
C      program "EULPLT".
C
      INCLUDE 'common_block_euler'
C
C
      DIMENSION  PSTAG(I_MAX,J_MAX), VMACH(I_MAX,J_MAX),
     &           PLOTVAR(I_MAX,J_MAX)
C
      OPEN(UNIT=7,FILE='euler.plt')
C
      WRITE(7,99) title
   99 FORMAT(A80)
C
      CP = RGAS*GAMMA/(GAMMA-1.)
      CV = CP/GAMMA
      WRITE(7,100) 1,NI,NJ,0,0,1,NI,0,0
      WRITE(7,101) CP,GAMMA

      DO 10 I=1,NI
      WRITE(7,101) (x(I,J),j=1,nj)
      WRITE(7,101) (y(i,j),j=1,nj)
  10  CONTINUE
  100 FORMAT(16I5)
  101 FORMAT(10F10.5)
C
C     CALCULATE THE SECONDARY VARIABLES
C
      DO I=1,NI
      DO J=1,NJ
      VX(I,J)  = ROVX(I,J)/RO(I,J)
      VY(I,J)  = ROVY(I,J)/RO(I,J)
      EKE      = 0.5*(VX(I,J)*VX(I,J) + VY(I,J)*VY(I,J))
      TSTAT    = (HSTAG(I,J) - EKE)/CP
      P(I,J)   = RO(I,J)*RGAS*TSTAT
      ROE(I,J) = RO(I,J)*(CV*TSTAT + EKE)
      END DO
      END DO
C
C    CALCULATE THE MACH NUMBER AND STAGNATION PRESSSURE
C
      DO I=1,NI
      DO J=1,NJ
      TSTAT = P(I,J)/RGAS/RO(I,J)
      VELSQ = VX(I,J)*VX(I,J) + VY(I,J)*VY(I,J)
      TSTAG = TSTAT + 0.5*VELSQ/CP
      VMACH(I,J) = SQRT(VELSQ/(GAMMA*RGAS*TSTAT))
      PSTAG(I,J) = P(I,J)*(TSTAG/TSTAT)**(1/FGA)
      END DO
      END DO
C
      WRITE(7,*) ' TIME STEP NUMBER', NSTEP
C
      WRITE(7,*)  ' AXIAL VELOCITY'
      DO 200 I=1,NI
 200  WRITE(7,201)(VX(I,J),J=1,NJ)
 201  FORMAT(10F10.4)
C
      WRITE(7,*)  ' Y  VELOCITY'
      DO 300 I=1,NI
 300  WRITE(7,201)(VY(I,J),J=1,NJ)
C
      WRITE(7,*)  ' RADIAL VELOCITY'
      DO 400 I=1,NI
 400  WRITE(7,201)(0.0,J=1,NJ)
C
      WRITE(7,*)  ' MACH NUMBER '
      DO 500 I=1,NI
 500  WRITE(7,201)(VMACH(I,J),J=1,NJ)
C
      WRITE(7,*)  ' STATIC PRESSURE'
      DO 600 I=1,NI
 600  WRITE(7,202)(P(I,J),J=1,NJ)
 202  FORMAT(10F10.1)
C
      WRITE(7,*)  ' DENSITY '
      DO 700 I=1,NI
 700  WRITE(7,201)(RO(I,J),J=1,NJ)
C
      WRITE(7,*)  ' VARIABLE PLOTVAR '
      DO 800 I=1,NI
 800  WRITE(7,201)(PLOTVAR(I,J),J=1,NJ)
C
      WRITE(7,*)  ' DEL ROVX '
      DO 900 I=1,NI
 900  WRITE(7,203)(DELROVX(I,J),J=1,NJ)
 203  FORMAT(10E10.4)
C
      CLOSE(7)
C
      RETURN
      END
