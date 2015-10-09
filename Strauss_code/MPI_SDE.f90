   program MPI_SDE
      
   
	include 'mpif.h'
        integer numtasks, rank, ierr, rc

!	-----------------------------------------------------------

	INTEGER :: i, j, k
	REAL :: E_begin, r_begin
	REAL :: j_BEGIN,counter,j_N, P,theta_begin


	INTEGER, PARAMETER :: Nk = 1
!	REAL :: Energies(Nk)


!	INTEGER :: file_counter = 100.0
	INTEGER :: task_counter = 0.0

	REAL, PARAMETER :: PI = 3.141592653589793

!	INTEGER :: SEEDS(7200)

!	REAL :: rgrid(400)


	REAL :: phi_begin
	REAL :: phi_earth(1)

    REAL :: r, th, ph
    REAL :: Br, Bph, B, B0


!	-----------------------------------------------------------
!	Different SEED for each process

!	CALL PARALLEL_SEEDS(SEEDS)
	
!	-----------------------------------------------------------

!   call MPI_INIT(ierr)
!   if (ierr .ne. MPI_SUCCESS) then
!      print *,'Error starting MPI program. Terminating.'
!      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
!   end if
!
!   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
!   print *, 'Number of processes:',numtasks,' Check process number . . .',rank

!	-----------------------------------------------------------

!	phi_begin = 0.0
!	E_begin = 0.004


!	Energies(1) = 0.001
!	Energies(2) = 0.002
!	Energies(3) = 0.005
!	Energies(4) = 0.01
!	Energies(5) = 0.02
!	Energies(6) = 0.05
!	Energies(7) = 0.1
!	Energies(8) = 0.2
!	Energies(9) = 0.5
!	Energies(10) = 1.0
!	Energies(11) = 2.0
!	Energies(12) = 5.0
!	Energies(13) = 10.0
!	Energies(14) = 20.0
!	Energies(15) = 50.0



	phi_earth(1) = 0.0



!	DO j = 2, 160, 1

!		phi_earth(j) = phi_earth(j - 1) + 2.0*PI/24.0



!	IF (phi_earth(j).GE.2.0*PI) THEN

 !     phi_earth(j) = phi_earth(j) - 2.0*PI

!	END IF

!Energies(1) = 0.10

!DO k = 2, NK, 1

!Energies(k) = Energies(k-1)*exp(0.1)

!END DO

!	END DO

!	rgrid(1) = 0.1

!	DO j = 2, 400, 1

!		rgrid(j) = rgrid(j - 1) + 0.34975

!	END DO


!rgrid(1) = 0.1
!rgrid(2) = 1.0
!rgrid(3) = 10.0
!rgrid(4) = 20.0
!rgrid(5) = 25.0
!rgrid(6) = 30.0
!rgrid(7) = 40.0
!rgrid(8) = 45.0
!rgrid(9) = 50.0
!rgrid(10) = 60.0
!rgrid(11) = 65.0
!rgrid(12) = 70.0
!rgrid(13) = 80.0
!rgrid(14) = 85.0
!rgrid(15) = 90.0


!	-----------------------------------------------------------

!        OPEN(700,file='Diffusion_Coefficients.txt',status='unknown')
!        OPEN(600,file='Solar_wind_density.txt',status='unknown')
!        OPEN(500,file='Alfven_speed.txt',status='unknown')
!        OPEN(400,file='Turbulence_parameters.txt',status='unknown')

!DO j = 1, 1, 1
!
!
!	DO i = 1, numtasks - 1, 1 
!
!		file_counter = file_counter + 1.0
!		task_counter = task_counter + 1
!
!		IF (rank.EQ.i) THEN
!
!WRITE(*,*) 'Started process number:', task_counter, 'on node no.:', rank
!
!	
!    E_begin = 0.100 !Energies(task_counter)
!
!    phi_begin = pi
!    theta_begin = pi/2.0
!    r_begin = 1.0

    B0 = 5.0/SQRT(2.0) * 0.06

    r = 1.0
    th = pi/2.0
    ph = 0.0

    CALL MAGNETICFIELD(r,th,ph,Br,Bph,B,B0)

    WRITE(*,*) 'Magnetic field at (r, th, ph) = (', r, ', ', th, ', ', ph, '): Br = ', Br, ', Bph = ', Bph,&
        ', B_Earth = ', B

    r = 10.0
    th = pi/3.0
    ph = pi/3.0

    CALL MAGNETICFIELD(r,th,ph,Br,Bph,B,B0)

    WRITE(*,*) 'Magnetic field at (r, th, ph) = (', r, ', ', th, ', ', ph, '): Br = ', Br, ', Bph = ', Bph

    r = 10.0
    th = -pi/3.0
    ph = 3.0*pi/2.0

    CALL MAGNETICFIELD(r,th,ph,Br,Bph,B,B0)

    WRITE(*,*) 'Magnetic field at (r, th, ph) = (', r, ', ', th, ', ', ph, '): Br = ', Br, ', Bph = ', Bph


!	Call different seed on each process
!			CALL sgrnd(SEEDS(i))
!
!	Call the code
    !CALL SDE(phi_begin,r_begin,E_begin,j_BEGIN,counter,j_N, P,theta_begin,SEEDS(task_counter))

!	Write the output to file
    !file_counter = 1.0
    !WRITE(file_counter,"(7(ES18.8))") r_begin, E_begin,j_BEGIN,counter, j_N*P*P, theta_begin, phi_begin

!	WRITE(*,*) 'Finished process number:', task_counter, 'with SEED:', SEEDS(task_counter)
    !WRITE(*,*) 'Output printed to file . . . ', file_counter


!		END IF
!
!	END DO ! the FOR
!
!END DO ! the FOR

!CLOSE(700)
!CLOSE(600)
!CLOSE(500)
!CLOSE(400)
!
!  call MPI_FINALIZE(ierr)

   end

!	-----------------------------------------------------------
!	-----------------------------------------------------------
!	-----------------------------------------------------------

 	SUBROUTINE SDE(phi_begin,r_begin,E_begin,j_BEGIN,counter,j_N, P,theta_begin,SEED)
		IMPLICIT NONE

!	-----------------------------------------------------------

	REAL :: RNumber1 = 0.0, RNumber2 = 0.0
	REAL :: RNumber3 = 0.0, RNumber4 = 0.0
	REAL, PARAMETER :: N = 10000

	DOUBLE PRECISION :: RANDOMNUMBERS
	INTEGER :: SEED

	REAL :: Counter
	REAL, PARAMETER :: r_inner = 0.01
	REAL :: E_begin, r_begin
	REAL, PARAMETER :: delta_T = 0.001
	REAL, PARAMETER :: r_boundary = 120.0, E_0 = 0.000511
!0.938
	
	
	REAL :: E_end, P
	REAL, PARAMETER :: PI = 3.141592653589793
	REAL :: phi_begin, phi, theta_begin, theta, theta_end

	REAl :: AA, CC, GG, MM

	REAL :: B11, B12, B13, B22, B23, B33
	
	REAL :: r, E, j_BEGIN, gamma, j_N

	INTEGER :: printer, endwhile

	REAL :: timespend

!	-----------------------------------------------------------
!       Files to write output

    OPEN(800,file='End_Position.txt',status='unknown')
	OPEN(999,file='End_energy.txt',status='unknown')
	OPEN(995,file='Time_spend.txt',status='unknown')
!	-----------------------------------------------------------


	RANDOMNUMBERS = 0.0

	Write(*,*) 'Energy: ', E_begin 
	WRITE(*,*) 'Radial distance: ', r_begin
	Write(*,*) 'Theta: ', theta_begin 
	WRITE(*,*) 'Phi: ', phi_begin


	printer = 0
	endwhile = 0

	Counter = 0.0
	j_N = 0.0
	j_BEGIN = 0.0


	DO WHILE (endwhile.EQ.0)

		r = r_begin
		E = E_begin
		phi = phi_begin
		theta = theta_begin

	timespend = 0.0
		
		DO WHILE ((r.LT.r_boundary).AND.(r.GT.r_inner))
		
				CALL RandomNumber(RNumber1,RNumber2,RANDOMNUMBERS,SEED)
				CALL RandomNumber(RNumber3,RNumber4,RANDOMNUMBERS,SEED)

			timespend = timespend + 1.0

	CALL Coefficients(r,phi,E,E_0,gamma,P,AA, CC, GG, delta_T, MM, theta, printer,B11, B12, B13, B22, B23, B33)
		
			r = r + AA*delta_T + B11*RNumber1 + B12*RNumber2 + B13*RNumber3

			theta = theta + MM*delta_T + B22*RNumber2 + B23*RNumber3
			
			phi = phi + CC*delta_T + B33*RNumber3

			E = E + GG*delta_T

!	Renormalize coordnates

	DO WHILE (theta.GT.PI)

      theta = 2.0*PI - theta
	  phi = phi - PI

	END DO

	DO WHILE (theta.LT.0.0)

      theta = -theta
	  phi = phi + PI

	END DO

	DO WHILE (phi.GT.2.0*PI)

		phi = phi - 2.0*PI

	END DO

	DO WHILE (phi.LT.0.0)

		phi = phi + 2.0*PI

	END DO


		END DO !End while loop  - particle trajectory


!	Check for LIS boundary


		IF (r.GT.r_begin) THEN

			E_end = E	
			Counter = Counter + 1.0
			theta_end = theta

	Write(*,*)  Counter

 		P = SQRT(E_end*(E_end + 2.0*E_0))

			CALL LIS(E_end,P,j_N,E_0)
			
			j_BEGIN = j_BEGIN + j_N

!	WRITE(800,"(3(ES18.8))") r, theta, phi
	WRITE(999,"(1(ES18.8))") E_end
	WRITE(995,"(1(ES18.8))") timespend*4.31*delta_T

		END IF

!		WRITE(*,*) counter

	IF ((Counter.GE.N)) THEN

		endwhile = 1

	END IF


	END DO !End while loop - number of particles



		P = SQRT(E_begin*(E_begin + 2.0*E_0))

		j_BEGIN = j_BEGIN/Counter*P*P



		CALL LIS(E_begin,P,j_N,E_0)

CLOSE(800)
CLOSE(999)
  CLOSE(995)
!	------------------------------------------------------

	RETURN
	END

!--------------------------------------------------------


!-----------------------------------------------------

	SUBROUTINE Coefficients(r,phi,E,E_0,gamma,P,AA, CC, GG, delta_T, MM, theta, printer,B11, B12, B13, B22, B23, B33)
		IMPLICIT NONE

	REAL :: r, E, E_0, V_sw, DivV_sw
	REAL :: BETA, gamma, P
	REAL :: phi, vdr, vdtheta, vdphi
	REAL :: delta_T, distance
	REAL, PARAMETER :: pi = 3.141592653589793

	REAL :: AA, CC, GG, MM, theta

	REAL :: B11, B12, B13, B22, B23, B33

!	Variables for diffusion tensor
	REAL :: K_parallel, lambda_parallel, K_perpr, K_perptheta
	REAL :: K_rr, K_rphi, K_phiphi, K_thetatheta, K_thetar, K_thetaphi
	REAL :: dK_rrdr, dK_rphi_dr, K_rthetadr
	REAL :: dK_thetardtheta, dKthetthet_dthet, dK_thetaphidtheta
	REAL :: dK_phirdphi, d_Kphiphidphi, dK_thetaphidphi

!	Variable for derivative calculations - r direction

	REAL :: delta_r = 0.01
	REAL :: r1
	REAL :: K_parallel1, lambda_parallel1, lambda_perpr1, K_perpr1, lambda_perptheta1, K_perptheta1
	REAL :: K_rr1, K_rphi1, K_phiphi1, K_thetatheta1, K_thetar1, K_thetaphi1

!	Variable for derivative calculations - theta direction

	REAL :: delta_theta = pi/180.0
	REAL :: theta2
	REAL :: K_parallel2, lambda_parallel2, lambda_perpr2, K_perpr2, lambda_perptheta2, K_perptheta2
	REAL :: K_rr2, K_rphi2, K_phiphi2, K_thetatheta2, K_thetar2, K_thetaphi2

!	Variable for derivative calculations - phi direction

	REAL ::	delta_phi = 2.0*pi/360.0
	REAL :: phi3
	REAL :: K_parallel3, lambda_parallel3, lambda_perpr3, K_perpr3, lambda_perptheta3, K_perptheta3
	REAL :: K_rr3, K_rphi3, K_phiphi3, K_thetatheta3, K_thetar3, K_thetaphi3

!	Variables for drift calculations
	REAL :: gamma_d, drift_coeff, pol_cycle = -1.0
        REAL :: charge = -1.0
	REAL :: Larmor, FS, Omega,sinpsi,cospsi

!	variables for the HMF field
	REAL :: Br, Bphi, B, Bearth
        REAL, PARAMETER :: B0 = 5.0/SQRT(2.0)*0.06

!       Variables for printing output
        INTEGER :: printer

!=======================================================================

!	Solar wind speed (400km/s)
	V_sw = 1.0

!	Divergence of solar wind (400km/s/AU)
	DivV_sw = 2.0*V_sw/r

	gamma = (E + 2.0*E_0)/(E + E_0)

	P = SQRT(E*(E + 2.0*E_0))
	BETA = P/SQRT(P*P + E_0*E_0)
	
!	----------------------------------------
	CALL MAGNETICFIELD(1.0,pi/2.0,0.0,Br,Bphi,B,B0)

	Bearth = B

!	Drift coefficients

!	larmor radius - AU
	CALL MAGNETICFIELD(r,theta,phi,Br,Bphi,B,B0)
	Larmor = 1.0/B*P/750.0

!	Drift reduction function
	FS = 0.5*10.0*P*P/(1.0 + 10.0*P*P)
!10.0*P*P/(1.0 + 10.0*P*P)

        CALL SPIRAL_ANGLE(r,theta,phi,Omega,sinpsi,cospsi)

	gamma_d = r*Omega*sin(theta)/V_sw
	drift_coeff = 2.0/3.0/pol_cycle/charge/B0*P*BETA*r/(1.0+gamma_d*gamma_d)/(1.0+gamma_d*gamma_d)*FS

	distance = abs(r*cos(theta))

	IF(theta.LT.pi/2.0)THEN

		vdr = drift_coeff*(-gamma_d/tan(theta))
		vdtheta = drift_coeff*(2.0*gamma_d + gamma_d*gamma_d*gamma_d)
		vdphi = drift_coeff*(gamma_d*gamma_d/tan(theta))

		ENDIF

	IF(theta.EQ.pi/2.0)THEN

		vdr = 0.0
		vdtheta = 0.0
		vdphi = 0.0

	ENDIF

	IF(theta.GT.pi/2.0)THEN
	
		vdr = -drift_coeff*(-gamma_d/tan(theta))
		vdtheta = -drift_coeff*(2.0*gamma_d + gamma_d*gamma_d*gamma_d)
		vdphi = -drift_coeff*(gamma_d*gamma_d/tan(theta))



	END IF
	

	IF(distance.LE.2.0*Larmor)THEN

    vdr = vdr + pol_cycle*charge*(0.457 - 0.412*distance/Larmor + 0.0915*distance*distance/Larmor/Larmor)*BETA*750.0*FS*sinpsi
    vdphi = vdphi + pol_cycle*charge*(0.457 - 0.412*distance/Larmor + 0.0915*distance*distance/Larmor/Larmor)*BETA*750.0*FS * cospsi

	ENDIF

!	----------------------------------------------
!	Diffusion coefficient

	CALL MAGNETICFIELD(r,theta,phi,Br,Bphi,B,B0)

    CALL DIFFUSION_COEFFICIENTS(Bearth,B,P,K_parallel,K_perpr,K_perptheta,BETA,r,theta,printer,V_sw)

!	Spherical coordinates	
	CALL SPHERICALCOORDINATES(r,theta,phi,K_parallel,K_perpr,K_perptheta,K_rr,&
K_rphi, K_phiphi, K_thetatheta, K_thetar, K_thetaphi)

!	--------------------------------------------
!	Derivatives to r

	r1 = r + delta_r
	
	CALL MAGNETICFIELD(r1,theta,phi,Br,Bphi,B,B0)

    CALL DIFFUSION_COEFFICIENTS(Bearth,B,P,K_parallel1,K_perpr1,K_perptheta1,BETA,r1,theta,printer,V_sw)

!	Spherical coordinates	
	CALL SPHERICALCOORDINATES(r1,theta,phi,K_parallel1,K_perpr1,K_perptheta1,K_rr1, &
K_rphi1, K_phiphi1, K_thetatheta1, K_thetar1, K_thetaphi1)

	dK_rrdr = (K_rr1 - K_rr)/delta_r
	dK_rphi_dr = (K_rphi1 - K_rphi)/delta_r
	K_rthetadr = (K_thetar1 - K_thetar)/delta_r

!	WRITE(*,*) dK_rrdr, K_parallel1

!	----------------------------------------------
!	Derivatives to theta

	theta2 = theta + delta_theta

	IF (theta2.GT.pi) THEN

		delta_theta = -1.0*delta_theta
		theta2 = theta + delta_theta

	ENDIF

	CALL MAGNETICFIELD(r,theta2,phi,Br,Bphi,B,B0)

        CALL DIFFUSION_COEFFICIENTS(Bearth,B,P,K_parallel2,K_perpr2,K_perptheta2,BETA,r,theta2,printer,V_sw)

!	Spherical coordinates	
	CALL SPHERICALCOORDINATES(r,theta2,phi,K_parallel2,K_perpr2,K_perptheta2,K_rr2, &
K_rphi2, K_phiphi2, K_thetatheta2, K_thetar2, K_thetaphi2)

	dK_thetardtheta = (K_thetar2 - K_thetar)/delta_theta
	dKthetthet_dthet = (K_thetatheta2 - K_thetatheta)/delta_theta
	dK_thetaphidtheta = (K_thetaphi2 - K_thetaphi)/delta_theta

!	----------------------------------------------
!	Derivatives to phi

	phi3 = phi + delta_phi
	
	IF (phi3.GT.2.0*pi) THEN

		delta_phi = -1.0*delta_phi
		phi3 = phi + delta_phi

	ENDIF
	
	CALL MAGNETICFIELD(r,theta,phi3,Br,Bphi,B,B0)

        CALL DIFFUSION_COEFFICIENTS(Bearth,B,P,K_parallel3,K_perpr3,K_perptheta3,BETA,r,theta,printer,V_sw)

!	Spherical coordinates	
	CALL SPHERICALCOORDINATES(r,theta,phi3,K_parallel3,K_perpr3,K_perptheta3,K_rr3, &
K_rphi3, K_phiphi3, K_thetatheta3, K_thetar3, K_thetaphi3)


	dK_phirdphi = (K_rphi3 - K_rphi)/delta_phi
	d_Kphiphidphi = (K_phiphi3 - K_phiphi)/delta_phi
	dK_thetaphidphi = (K_thetaphi3 - K_thetaphi)/delta_phi

!	----------------------------------------------
!	Numerical coefficients

!	First order derivatives

	GG = DivV_sw/3.0*gamma*E

	AA = 2.0/r*K_rr + dK_rrdr - V_sw - vdr + 1.0/r/tan(theta)*K_thetar + 1.0/r*dK_thetardtheta &
			+ 1.0/r/sin(theta)*dK_phirdphi
	
	MM = 1.0/r/r*K_thetar + 1.0/r*K_rthetadr + 1.0/r/r/tan(theta)*K_thetatheta &
			+ 1.0/r/r*dKthetthet_dthet + 1.0/r/r/sin(theta)*dK_thetaphidphi - vdtheta/r

	CC = K_rphi/r/r/sin(theta) + dK_rphi_dr/r/sin(theta) + dK_thetaphidtheta/r/r/sin(theta) &
			+ d_Kphiphidphi/r/r/sin(theta)/sin(theta) - vdphi/r/sin(theta)

!	Tensor elements

	B11 = SQRT(2.0*delta_T)*SQRT((K_phiphi*K_thetar*K_thetar - 2.0*K_rphi*K_thetar*K_thetaphi + K_rr*K_thetaphi*K_thetaphi + &
					K_thetatheta*K_rphi*K_rphi - K_rr*K_phiphi*K_thetatheta)/(K_thetaphi*K_thetaphi - K_phiphi*K_thetatheta))

	B12 = SQRT(2.0*delta_T)*((K_rphi*K_thetaphi - K_thetar*K_phiphi)/(K_thetaphi*K_thetaphi - K_thetatheta*K_phiphi)* &
					SQRT(K_thetatheta - K_thetaphi*K_thetaphi/K_phiphi))

	B13 = SQRT(2.0*delta_T)*(K_rphi/SQRT(K_phiphi))

	B22 = SQRT(2.0*delta_T)*(SQRT(K_thetatheta - K_thetaphi*K_thetaphi/K_phiphi)/r)

	B23 = SQRT(2.0*delta_T)*(K_thetaphi/r/SQRT(K_phiphi))

	B33 = SQRT(2.0*delta_T)*(SQRT(K_phiphi)/r/sin(theta))

	RETURN
	END

!	=================================================
SUBROUTINE MAGNETICFIELD(r,theta,phi,Br,Bphi,B,B0)
	IMPLICIT NONE

	REAL, PARAMETER :: pi = 3.141592653589793
	REAL :: V_sw = 1.0
	REAL :: R_Sun = 0.005, Omega1 = 2.0*pi/25.4/3600.0/24.0, PROTIME = 1.496d8/400.0
	REAL :: Omega = 0.0, B0
	REAL :: Br, Bphi,B
	REAL :: r, theta, phi

	Omega = Omega1*PROTIME

	Br = B0/r/r
	Bphi = -B0/r/r*Omega*(r - R_Sun)*sin(theta)/V_sw

	B = SQRT(Br*Br + Bphi*Bphi)

RETURN
END
!	=================================================
SUBROUTINE SPIRAL_ANGLE(r,theta,phi,Omega,sinpsi,cospsi)
	IMPLICIT NONE

    REAL :: r, theta, phi, V_sw = 1.0
    REAL :: Omega, sinpsi, cospsi, tanpsi
    REAL, PARAMETER :: pi = 3.141592653589793
    REAL :: R_Sun = 0.005, Omega1 = 2.0*pi/25.4/3600.0/24.0, PROTIME = 1.496d8/400.0

    Omega = Omega1*PROTIME
    tanpsi = Omega*(r - R_Sun)*sin(theta)/V_sw
	cospsi = 1.0/SQRT(1.0 + tanpsi*tanpsi)
	sinpsi = SQRT(1.0 - cospsi*cospsi)

RETURN
END

!==========================================================
SUBROUTINE DIFFUSION_COEFFICIENTS(Bearth,B,P,K_parallel,K_perpr,K_perptheta,BETA,r,theta,printer,V_sw)

       REAL, PARAMETER :: pi = 3.141592653589793

       REAL :: Bearth, B, P, BETA, V_sw
       REAL :: K_parallel,K_perpr,K_perptheta
       REAL :: lambda_parallel, r, theta, lambda_perp
       INTEGER :: printer

       REAL :: density, alfven, kmin, l2D
       REAL :: larmor, delta_slab, delta_2D, q, Rc

       REAL :: phi, Omega, sinpsi, cospsi, pin
       REAL :: Qc, bc, aa, bb, kii, kd, OmegaCi


       CALL SPIRAL_ANGLE(r,theta,phi,Omega,sinpsi,cospsi)

!      ----- Solar wind mass density
       density = 5.0d6/r/r/V_sw*1.67d-27

!      ----- Alfven speed
       alfven = B/0.06*1.0d-9/SQRT(density*1.26d-6)/1000.0/400.0

!      ----- kmin - break energy inertial range
       kmin = 32.0/SQRT(r)

!      ----- larmor radius
       larmor = 1.0/B*P/750.0

!      ----- slab variance
       delta_slab = 13.2*r**(-2.7)*0.06*0.06

!      ----- 2D variance
       delta_2D = 4.0*delta_slab

!      ----- Kolmagorov index
       q = 5.0/3.0
!      ----- Dissipation range index
       pin = 2.6

       Rc = kmin*larmor

!      ---- 2D correlation length
       l2D = 3.1d-3*SQRT(r)

!      ----- proton gyrofrequency
       OmegaCi = 1.6d-19*B/0.06*1.0d-9/1.672d-27

!      ----- ion inertial scale
       kii = 2.0*pi*OmegaCi*sinpsi/alfven

!      ----- kd - break inertial dissipation range
       aa = 0.152
       bb = 0.451
       kd = 2.0*pi/V_sw*(aa + bb/2.0/pi*kii*V_sw)*374000.0

       Qc = kd*larmor
       bc = BETA*750.0/2.0/alfven

!     -----------------------------------------------
!     lambda parallel
!     protons
      lambda_parallel = 3.0*q/pi*Rc*Rc/kmin*B*B/delta_slab*(0.25 + 2.0/(2.0-q)/(4.0-q)/Rc**(q))

!     electrons
!      lambda_parallel = 3.0*q/SQRT(pi)/(q-1.0)*Rc*Rc/bc/kmin*B*B/delta_slab*( &
!                                                                             bc/4.0/SQRT(pi)+&
!                                                                             2.0/SQRT(pi)/(2.0-q)/(4.0-q)*bc/Rc**(q)+&
!                                                                             (1.0/0.886+1.0/SQRT(pi)/1.5)*bc**(pin-1.0)/Rc**(q)/Qc**(pin-q))




!      lambda_parallel = 0.05*Bearth/B*P

!    IF (P.LT.1.0) THEN

!    lambda_parallel = 0.05*Bearth/B

!    END IF

!    K_parallel = lambda_parallel*BETA*750.0/3.0


    K_parallel = 25.0*P*(1.0 + r)

    IF (P.LT.1.0) THEN

    K_parallel = 25.0*1.0*(1.0 + r)

    END IF

!    -----------------------------------------------
!    lambda perp

!     lambda_perp = (1.0/SQRT(3.0)/2.0/SQRT(pi)*l2D*1.128/2.565*delta_2D/B/B)**(2.0/3.0)*lambda_parallel**(1.0/3.0)*r**(0.65)


    K_perpr = 0.02*K_parallel
! lambda_perp*BETA*750.0/3.0
	K_perptheta = K_perpr

!	----------------------------------------------
!	Write some output

	IF (printer.EQ.0) THEN

		WRITE(700,"(5(ES18.8))") r, P, theta, lambda_parallel, lambda_perp
		WRITE(600,"(4(ES18.8))") r, theta, density, density/1.0d6/1.67d-27
		WRITE(500,"(4(ES18.8))") r, theta, alfven*400.0, B/0.06
		WRITE(400,"(10(ES18.8))") r, theta, kmin, larmor, delta_slab, delta_2D, l2D, kd, kii, OmegaCi

		printer = printer + 1

	END IF

RETURN
END
!-------------------------------------------------------
	SUBROUTINE SPHERICALCOORDINATES(r,theta,phi,K_parallel,K_perpr,K_perptheta,K_rr, &
K_rphi, K_phiphi, K_thetatheta, K_thetar, K_thetaphi)
		IMPLICIT NONE

	REAL, PARAMETER :: pi = 3.141592653589793
	REAL, PARAMETER :: V_sw = 1.0

		REAL :: r,theta,phi
		REAL :: K_parallel,K_perpr,K_perptheta
		REAL :: K_rr, K_rphi, K_phiphi, K_thetatheta, K_thetar, K_thetaphi

	REAL :: tanpsi, cospsi, sinpsi
	REAL :: R_Sun = 0.005, Omega1 = 2.0*pi/25.4/3600.0/24.0, PROTIME = 1.496d8/400.0
	REAL :: Omega = 0.0

	Omega = Omega1*PROTIME

	tanpsi = Omega*(r - R_Sun)*sin(theta)/V_sw
	cospsi = 1.0/SQRT(1.0 + tanpsi*tanpsi)
	sinpsi = SQRT(1.0 - cospsi*cospsi)

			K_rr = K_parallel*cospsi*cospsi + K_perpr*sinpsi*sinpsi
			K_phiphi = K_parallel*sinpsi*sinpsi + K_perpr*cospsi*cospsi
			K_rphi = (K_perpr - K_parallel)*cospsi*sinpsi
			K_thetatheta = K_perptheta

	!	For the parker HMF
	K_thetar = 0.0
	K_thetaphi = 0.0



	RETURN
	END
!##############################################################

subroutine anaptensor(r,t,p,kp,ks1,ks2,k11,k13,k33,k22,k12,k23)

!	Called in this order!!
!	(r,theta,phi,K_parallel,K_perpr,K_perptheta,K_rr, K_rphi, K_phiphi, K_thetatheta, K_thetar, K_thetaphi)

!	subroutine anaptensor(r,t,p,ks1,ks2,kp,k11,k12,k13,k22,k23,k33)

  implicit none

!	Changed all double precision to real kind

  REAL :: r,p,t,ks1,ks2,kp,dr,dp,dt

  REAL :: k11,k12,k13,k22,k23,k33

  REAL :: DD,EE,FFF,rs,tanx,tanr,tant,tanp

  REAL :: nanr,nant,nanp,banr,bant,banp

!  parameter(pi=3.1415926535897932384626433832795028841972)
  REAL, PARAMETER :: pi = 3.141592653589793

!	DD = omega/V_sw (in units of AU)

  !analytic comparison
  DD=1.06964
  rs=0.005

  tanx=DD*(r-rs)*sin(t)
  
  EE=tanx/r+DD*sin(t)/(1.+tanx**2.)
  FFF=tanx**2.*cos(t)/(r*sin(t))

  tanr=1./sqrt(1.+tanx**2.)
  tant=0.0
  tanp=-tanx/sqrt(1.+tanx**2.)
  
  nanr=-EE*tanx/sqrt(EE**2.*(tanx**2.+1.)+FFF**2.)
  nant=-FFF/sqrt(EE**2.*(tanx**2.+1.)+FFF**2.)
  nanp=-EE/sqrt(EE**2.*(tanx**2.+1.)+FFF**2.)
  
  !write(*,*) sqrt(nanr**2.+nant**2.+nanp**2.)
  !write(*,*) sqrt(nrn**2.+ntn**2.+npn**2.)
  
  banr=-FFF*tanx/sqrt(FFF**2.*(tanx**2.+1.) + EE**2.*(tanx**2.+1.)**2.)
  bant=EE*(1.+tanx**2.)/sqrt(FFF**2.*(tanx**2.+1.) + EE**2.*(tanx**2.+1.)**2.)
  banp=-FFF/sqrt(FFF**2.*(tanx**2.+1.) + EE**2.*(tanx**2.+1.)**2.)

  k11=ks1*nanr**2. + ks2*banr**2. + kp*tanr**2.
  k12=ks1*nanr*nant + ks2*banr*bant + kp*tanr*tant
  k13=ks1*nanr*nanp + ks2*banr*banp + kp*tanr*tanp
  k22=ks1*nant**2. + ks2*bant**2 + kp*tant**2
  k23=ks1*nant*nanp + ks2*bant*banp + kp*tant*tanp
  k33=ks1*nanp**2. + ks2*banp**2. + kp*tanp**2.

  return
end
!-----------------------------------------------------
	SUBROUTINE LIS(E_end,P,j_N,E_0)
		IMPLICIT NONE

	REAL :: E_end,P,j_N
	REAL :: j_LIS
	REAL :: theta, E_0, a, b, c, d


		P = SQRT(E_end*(E_end + 2.0*E_0))

!	Protons
!		IF (E_end.LT.1.0) THEN
 
!			j_LIS = exp(4.64 - 0.08*log(E_end)*log(E_end) &
!				-2.91*SQRT(E_end))/P/P

!		ELSE
	
!			j_LIS = exp(3.22 - 2.86*log(E_end) - 1.5/E_end)/P/P

!		END IF

!		j_N = j_LIS

!        Electron
		IF (P.LT.0.0026) THEN

				a=126.067948848145
				b=0.2567973348983205
				c=1.95129069032843
				d=0.0171199701826333
				
				j_LIS = (a + c*log(P))/(1.0 + b*log(P) + d*log(P)**2)*1.7/P/P
 
!			j_LIS = (214.32+3.32*log(P))/(1.0+0.26*log(P)+0.02*log(P)*log(P))/P/P

		END IF

		IF ((P.GE.0.0026).AND.(P.LT.0.1)) THEN
 
			j_LIS = 1.7*((52.55+23.01*P)/(1.0+148.62*P))**(2)/P/P

		END IF

		IF ((P.GE.0.1).AND.(P.LE.10.0)) THEN
 
			j_LIS = (1555.89+17.36*P-3.4d-3*P*P+5.13d-7*P*P*P)/(1.0-11.22*P+7532.93*P*P+2405.01*P*P*P+103.87*P*P*P*P)/P/P

		END IF

		IF (P.GT.10.0) THEN
 
			j_LIS = 1.7*exp(-0.89-3.22*log(P))/P/P

		END IF



		j_N = j_LIS



	RETURN
	END
!-----------------------------------------------------
	SUBROUTINE RandomNumber(RNumber1,RNumber2,RANDOMNUMBERS,SEED)
	
	INTEGER :: j
	INTEGER, PARAMETER :: no = 1000000
	DOUBLE PRECISION :: Random1(0:7), Random2(0:7)
	DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793
	REAL :: RNumber1, RNumber2
	DOUBLE PRECISION :: RANDOMNUMBERS, GRND
	INTEGER :: SEED
	
	Random1(mod(0,8))=grnd(SEED)
        Random2(mod(0,8))=grnd(SEED)

DO WHILE (Random1(0).EQ.0.0)

Random1(mod(0,8))=grnd(SEED) 

END DO

DO WHILE (Random2(0).EQ.0.0)

Random2(mod(0,8))=grnd(SEED) 

END DO

	RNumber1 = SQRT(-2.0*LOG(Random1(0)))*COS(2.0*PI*Random2(0))
	RNumber2 = SQRT(-2.0*LOG(Random1(0)))*SIN(2.0*PI*Random2(0))
	
	RANDOMNUMBERS = RANDOMNUMBERS + 2.0

	RETURN
	END
!-----------------------------------------------------
      subroutine sgrnd(seed)
!
      implicit integer(a-z)
!
! Period parameters
      parameter(N     =  624)
!
      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)


!		WRITE(*,*) SEED
!		WRITE(*,*) 'mt(mti)', mt(mti)



 1000 continue
!
      return
      end
!----------------------------------------------------
      double precision function grnd(SEED)
!
      implicit integer(a-z)
!
! Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
!                                    constant vector a
      parameter(UMASK = -2147483647)
!                                    most significant w-r bits
      parameter(LMASK =  2147483647)
!                                    least significant r bits
! Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
!
      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
!                     mti==N+1 means mt[N] is not initialized
!
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!                        mag01(x) = x ! MATA for x=0,1
!
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
!
      if(mti.ge.N) then
!                       generate N words at one time
        if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
          call sgrnd(SEED)
!                              a default initial seed is used
        endif
!
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
!
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
!
      return
      end

!===================================================================
	SUBROUTINE PARALLEL_SEEDS(SEEDS)

	IMPLICIT NONE

	INTEGER :: SEEDS(7200)

	INTEGER :: n
	INTEGER :: k

	OPEN(unit = 2, file = "Seedlings.txt")

	DO k = 1, 7200, 1

	READ(2,*) n

	SEEDS(k) = n

	END DO

	WRITE(*,*) 'Parallel seeds initiated . . .'

	CLOSE(2)

	RETURN
	END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
