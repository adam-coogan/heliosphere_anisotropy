 	PROGRAM SDE
    USE IFPORT
    IMPLICIT NONE

!	-----------------------------------------------------------
!	Initiate variables
!	SEED_0 - SEED for random number generator
!	N = number of iterations

	DOUBLE PRECISION :: SEED_1 = 975635
	DOUBLE PRECISION :: RNumber1 = 0.0, RNumber2 = 0.0
    CHARACTER(len=32) :: Nstr
	INTEGER :: N
    CHARACTER(len=32) :: runFileName

	DOUBLE PRECISION :: RANDOMNUMBERS

	INTEGER, PARAMETER :: Nk = 1, Nr = 1
!100
	REAL :: Energies(Nk), radial_grid(Nr)
	INTEGER :: i,k
	REAL :: Counter, delta_E = 0.3
	REAL, PARAMETER :: r_inner = 0.005
    CHARACTER(len=32) :: E_beginstr
	REAL :: E_begin, r_begin
	REAL, PARAMETER :: delta_T = 0.001
	REAL, PARAMETER :: r_boundary = 140.0, E_0 = 0.000511
	REAL :: E_end, P
	REAL, PARAMETER :: PI = 3.141592653589793
	REAL :: phi_begin = 0, phi, theta_begin = PI/2.0, theta, theta_end

	REAl :: AA, BB, CC, DD, GG, HH, MM, NN
	
	REAL :: r, E, j_BEGIN, gamma, j_N
    INTEGER NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!	Jovian stuff

	REAL, PARAMETER :: PROTIME = 1.496d8/400.0
	REAL, PARAMETER :: Omega_jupiter = 2.0*pi/4333.0/3600.0/24.0*PROTIME
	REAL :: phi_jupiter, theta_jupiter = pi/2.0
	REAL, PARAMETER :: r_jupiter = 5.2
	REAL, PARAMETER :: delta_phi_jupiter = 0.009*2.0
!0.0046
	REAL, PARAMETER :: delta_theta_jupiter = 0.009*2.0
!0.0046
	REAL :: rbegin_jupiter = r_jupiter - 0.0477*2.0
	REAL :: rend_jupiter = r_jupiter + 0.095*2.0
	REAL :: counter_jupiter, j_JUP, j_JUP2
	INTEGER :: jovian_source

	INTEGER :: printer, endwhile

	REAL :: timespend

    ! Get command line arguments:
    !   Number of runs
    !   Energy at Earth
    !   Output filename relative to rundata, without .csv extension.
    CALL get_command_argument(1, E_beginstr)
    READ(E_beginstr, *) E_begin
    CALL get_command_argument(2, Nstr)
    READ(Nstr, *) N
    CALL get_command_argument(3, runFileName)

    WRITE(*,*) "N = ", N
    WRITE(*,*) "E_begin = ", E_begin, " GeV"
    WRITE(*,*) 'rundata/' // TRIM(runFileName) // '.csv'

!	Output files

	!OPEN(100,file='Standard_Deviates.txt',status='unknown')
	!OPEN(200,file='ParticleTrajectory1.txt',status='unknown')
	!OPEN(300,file='ParticleTrajectory2.txt',status='unknown')	
	!OPEN(400,file='ParticleTrajectory3.txt',status='unknown')
	!OPEN(500,file='DifferentialIntensity.txt',status='unknown')
	!OPEN(600,file='PlasmaProfiles.txt',status='unknown')
	!OPEN(700,file='DiffusionCoefficients.txt',status='unknown')
	OPEN(800,file='rundata/' // TRIM(runFileName) // '.csv',status='unknown')
	!Open(900,file='Jovian-phi.txt',status='unknown')
	!OPEN(999,file='SEEDS.txt',status='unknown')
!	-----------------------------------------------------------

! Seed the random number generator
!CALL RANDOM_SEED()
! Generates new seed for each invocation of the program
CALL init_random_seed()

!	Energies(1) = 0.001
!   Energies(2) = 0.002
!   Energies(3) = 0.005

!	Energies(1) = 0.100
!	Energies(2) = 0.02
!	Energies(3) = 0.05
!	Energies(4) = 0.1
!	Energies(5) = 0.2
!	Energies(6) = 0.3
!	Energies(7) = 0.5
!	Energies(8) = 1.0
!	Energies(9) = 2.0
!	Energies(12) = 5.0
!	Energies(13) = 10.0

RANDOMNUMBERS = 0.0


!DO k = 2, Nk, 1

!Energies(k) = Energies(k - 1)*exp(delta_E)

!END DO

radial_grid(1) = 1.0

!	DO i = 2, Nr, 1

!		radial_grid(i) = radial_grid(i - 1) + 0.3

!	END DO
	
r_begin = 1.0d0
theta_begin = PI / 2.0d0
phi_begin = 0.0d0

WRITE(*,*) 'Energy: ', E_begin 

printer = 0
endwhile = 0

Counter = 0.0
Counter_jupiter = 0.0
j_N = 0.0
j_BEGIN = 0.0
j_JUP2 = 0.0
j_JUP = 0.0

DO WHILE (endwhile.EQ.0)
    !	Write(*,*)  Counter_jupiter, ' ::', Counter, ' ::',N, k, '  ::', Nk
    r = r_begin
    E = E_begin
    phi = phi_begin
    theta = theta_begin
    
    jovian_source = 0
    phi_jupiter = PI
    
    timespend = 0.0

    DO WHILE ((r.LT.r_boundary).AND.(r.GT.r_inner).AND.(jovian_source.EQ.0)) ! Particle trajectory
        CALL RandomNumber(RNumber1,RNumber2,SEED_1,RANDOMNUMBERS)

        ! TODO: DEBUG!!!
!        WRITE(*,*) "s = ", timespend*delta_T*4.3287*86400

        CALL Coefficients(r,phi,E,E_0,gamma,P,AA, BB, CC, DD, GG, HH, delta_T, MM, NN, theta, printer)

        r = r + AA*delta_T + BB*RNumber1 + HH*RNumber2
        
        phi = phi + CC*delta_T + DD*RNumber2
        
        CALL RandomNumber(RNumber1,RNumber2,SEED_1,RANDOMNUMBERS)
        
        theta = theta + MM*delta_T + NN*RNumber2
        
        E = E + GG*delta_T

        ! TODO DEBUG!!!
!        WRITE(*,*) "dr_ds = ", AA
!        WRITE(*,*) "dr_dWr = ", BB / sqrt(delta_T)
!        WRITE(*,*) "dr_dWph = ", HH / sqrt(delta_T)
!        WRITE(*,*) "dth_ds = ", MM
!        WRITE(*,*) "dth_dWth = ", NN / sqrt(delta_T)
!        WRITE(*,*) "dph_ds = ", CC
!        WRITE(*,*) "dph_dWph = ", DD / sqrt(delta_T)
!        WRITE(*,*) "dE_ds = ", GG
!        WRITE(*,*) ""

        phi_jupiter = phi_jupiter - Omega_jupiter*delta_T
        
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
        
        !	Check Jupiter position
        IF ((r.GT.rbegin_jupiter).AND.(r.LT.rend_jupiter)) THEN
            IF ((phi.GT.(phi_jupiter - delta_phi_jupiter)).AND.(phi.LT.(phi_jupiter + delta_phi_jupiter))) THEN
            IF ((theta.GT.(theta_jupiter - delta_theta_jupiter)).AND.(theta.LT.(theta_jupiter + delta_theta_jupiter))) THEN
                    Counter_jupiter = Counter_jupiter + 1.0
                    jovian_source = 1
                    E_end = E
                    
                    CALL JOVIAN_LIS(E_end, E_0, j_JUP)
                    
                    j_JUP2 = j_JUP2 + j_JUP
                    
                    !WRITE(*,*) Counter, Counter_jupiter, N
                    !WRITE(*,*) k, ':', Nk 
                    
!                            WRITE(900,"(2(ES18.8))") abs(phi_jupiter - PI)*180.0/PI, timespend*delta_T*4.3287
!                            WRITE(999,"(1(ES18.8))") RANDOMNUMBERS
                END IF
            END IF
        END IF
        
        timespend = timespend + 1.0
    END DO !End while loop  - particle trajectory
    
    !	Check for LIS boundary
    IF ((r.GT.r_begin).AND.(jovian_source.EQ.0)) THEN
        Counter = Counter + 1.0
        
        !!!!! Changed to my output format!
        !WRITE(800,"(7(ES18.8))") counter, r, E, theta, phi, E_end, timespend*delta_T*4.3287
        WRITE(800,"(5(ES18.8))") r, theta, phi, E, timespend*delta_T*4.3287*86400
    END IF

    IF (MOD(Counter, 100.0).EQ.0.0) THEN
        WRITE(*,*) "Percent complete = ", Counter/N * 100.0
    END IF
    
    IF ((Counter.GT.N)) THEN
        endwhile = 1
    END IF

END DO !End while loop - number of particles

P = SQRT(E_begin*(E_begin + 2.0*E_0))

WRITE(*,*) "j_BEGIN = ", j_BEGIN
j_BEGIN = j_BEGIN/Counter*P*P
j_JUP2 = j_JUP2/Counter_jupiter*P*P

CALL LIS(E_begin,P,j_N,E_0)

WRITE(500,"(8(ES18.8))") r_begin, E_begin,j_BEGIN,counter, j_N*P*P, theta_begin, j_JUP2, Counter_jupiter


	Close(100)
	Close(200)
	Close(300)
	Close(400)
	Close(500)
	Close(600)
	Close(700)
	Close(800)
	Close(900)
	Close(999)

!	------------------------------------------------------

	END PROGRAM SDE

!--------------------------------------------------------
	SUBROUTINE RandomNumber(RNumber1,RNumber2,SEED_1,RANDOMNUMBERS)
		IMPLICIT NONE

	REAL :: SDeviate1, SDeviate2
	DOUBLE PRECISION :: RNumber1, RNumber2
	DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793

	DOUBLE PRECISION :: a = 7.0**(5)
	DOUBLE PRECISION :: M = 2.0**(101) - 1.0
	DOUBLE PRECISION :: c = 0.0
	DOUBLE PRECISION :: Random1, Random2
	DOUBLE PRECISION :: SEED_1
	DOUBLE PRECISION :: RANDOMNUMBERS

! Strauss' really old-school method
!	Random1 = MOD(a*SEED_1 + c,M)
!	Random2 = MOD(a*Random1 + c,M)
!	SEED_1 = Random2
!	SDeviate1 = Random1/M
!	SDeviate2 = Random2/M
!	RNumber1 = SQRT(-2.0*ALOG(SDeviate1))*COS(2.0*PI*(SDeviate2))
!	RNumber2 = SQRT(-2.0*ALOG(SDeviate1))*SIN(2.0*PI*(SDeviate2))

! The smooth way of doing it
    CALL RANDOM_NUMBER(SDeviate1)
    CALL RANDOM_NUMBER(SDeviate2)
	RNumber1 = SQRT(-2.0*ALOG(1.0 - SDeviate1))*COS(2.0*PI*(1.0 - SDeviate2))
	RNumber2 = SQRT(-2.0*ALOG(1.0 - SDeviate1))*SIN(2.0*PI*(1.0 - SDeviate2))

    ! TODO: DEBUG!!!
!    WRITE(*,*) "seed = ", SEED_1
!    WRITE(*,*) "r1, r2 = ", RNumber1, RNumber2

	RANDOMNUMBERS = RANDOMNUMBERS + 2.0
   !WRITE(999,"(1(ES18.8))") SEED_1

	RETURN
	END

    ! Generates a different seed each time the program is run
    subroutine init_random_seed()
        use iso_fortran_env, only: int64
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
        INTEGER*4 :: getpid

        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(t)
            if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24_int64 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            end if
            pid = GETPID()
            t = ieor(t, int(pid, kind(t)))
            do i = 1, n
                seed(i) = lcg(t)
            end do
        end if
        call random_seed(put=seed)
    contains
        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)
            integer :: lcg
            integer(int64) :: s
            if (s == 0) then
            s = 104729
            else
            s = mod(s, 4294967296_int64)
            end if
            s = mod(s * 279470273_int64, 4294967291_int64)
            lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_random_seed

!-----------------------------------------------------

	SUBROUTINE Coefficients(r,phi,E,E_0,gamma,P,AA, BB, CC, DD, GG, HH, delta_T, MM, NN, theta, printer)
		IMPLICIT NONE

	REAL :: r, E, E_0, V_sw, DivV_sw, K_parallel
	REAL :: lambda_parallel, BETA, gamma, P
	REAL :: phi, vdr, vdtheta, vdphi
	REAL :: lambda_perpr, K_perpr, K_rr, K_rphi, dK_rrdr, dK_rphi_dr, K_phiphi
	REAL :: delta_T, distance, dKthetthet_dthet
	REAL, PARAMETER :: pi = 3.141592653589793
	INTEGER :: printer

	REAL :: AA, BB, CC, DD, GG, HH

	REAL :: tanpsi, cospsi, sinpsi
	REAL :: R_Sun = 0.005, Omega1 = 2.0*pi/25.4/3600.0/24.0, PROTIME = 1.496d8/400.0
	REAL :: Omega = 0.0
	REAL :: MM, NN, theta

	REAL :: lambda_perptheta, K_perptheta, K_thetatheta, F_theta

!	Variable for gradient calculations

	REAL, PARAMETER :: delta_r = 0.01
	REAL :: r2, K_parallel2, K_perpr2, lambda_parallel2, lambda_perpr2
	REAL :: K_rr2, K_rphi2, tanpsi2, cospsi2, sinpsi2

!	Variables for drift calculations

	REAL :: gamma_d, drift_coeff, pol_cycle = 1.0
	REAL, PARAMETER :: B0 = 5.0*0.06 ! Due to weird B-field units
	REAL :: Bmag, Larmor, Krphidphi, FS

	Omega = Omega1*PROTIME

!	Solar wind speed (400km/s)
	V_sw = 1.0


!	Divergence of solar wind (400km/s/AU)
	DivV_sw = 2.0*V_sw/r

	gamma = (E + 2.0*E_0)/(E + E_0)

	P = SQRT(E*(E + 2.0*E_0))
	BETA = P/SQRT(P*P + E_0*E_0)
	
!	GG = coefficnect for the energy term
	GG = DivV_sw/3.0*gamma*E

!	magnetic field magnitude
	tanpsi = Omega*(r - R_Sun)*sin(theta)/V_sw
	cospsi = 1.0/SQRT(1.0 + tanpsi*tanpsi)
	sinpsi = SQRT(1.0 - cospsi*cospsi)
	
	Bmag = B0/SQRT(1.0 + (1.0 - R_sun)*(1.0 - R_Sun))/r/r/cospsi

!	larmor radius - AU
	Larmor = 1.0/Bmag*P/750.0

!	Drift reduction function
	FS = 10.0*P*P/(1.0 + 10.0*P*P)

!	Drift coefficients
	gamma_d = r*Omega*sin(theta)/V_sw
	drift_coeff = 2.0/3.0/pol_cycle/B0*P*BETA*r/(1.0+gamma_d*gamma_d)/(1.0+gamma_d*gamma_d)*FS

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

		vdr = vdr + pol_cycle*(0.457 - 0.412*distance/Larmor + 0.0915*distance*distance/Larmor/Larmor)*BETA*750.0*FS

	ENDIF

!	----------------------------------------------
!	Diffusion coefficienst
	IF (P.LT.1.0) THEN

		lambda_parallel = 0.15*(1.0 + r)

	ELSE

		lambda_parallel = 0.15*(1.0 + r)*P

	END IF

!	Parallel
	K_parallel = lambda_parallel*BETA*750.0/3.0

!	Perp = r
	K_perpr = 0.01*K_parallel
	lambda_perpr = K_perpr*3.0/beta/750.0

!	Perp - theta
	K_perptheta = 0.01*K_parallel

!	Spherical coordinates
    K_rr = K_parallel*cospsi*cospsi + K_perpr*sinpsi*sinpsi
    K_phiphi = K_parallel*sinpsi*sinpsi + K_perpr*cospsi*cospsi
    K_rphi = (K_perpr - K_parallel)*cospsi*sinpsi
    K_thetatheta = K_perptheta

!	Parker HMF - Symmetry about phi
	Krphidphi = 0.0

!	Numerical coefficients
	BB = SQRT(2.0*K_rr*delta_T - 2.0*delta_T*K_rphi*K_rphi/K_phiphi)
	DD = SQRT(2.0*K_phiphi*delta_T)/r/sin(theta)
	HH = sqrt(2.0)*K_rphi/SQRT(K_phiphi)*sqrt(delta_T)

!	----------------------------------------------
!	r - Derivatives of the diffusion coefficients

	r2 = r + delta_r

			tanpsi2 = Omega*(r2 - R_Sun)*sin(theta)/V_sw
			cospsi2 = 1.0/SQRT(1.0 + tanpsi2*tanpsi2)
			sinpsi2 = SQRT(1.0 - cospsi2*cospsi2)


	IF (P.LT.1.0) THEN

		lambda_parallel2 = 0.15*(1.0 + r2)

	ELSE

		lambda_parallel2 = 0.15*(1.0 + r2)*P

	END IF

!	Parallel
	K_parallel2 = lambda_parallel2*BETA*750.0/3.0

!	Perp = r
	K_perpr2 = 0.01*K_parallel2



	K_rr2 = K_parallel2*cospsi2*cospsi2 + K_perpr2*sinpsi2*sinpsi2
	K_rphi2 = (K_perpr2 - K_parallel2)*cospsi2*sinpsi2

	dK_rrdr = (K_rr2 - K_rr)/delta_r
	dK_rphi_dr = (K_rphi2 - K_rphi)/delta_r
		
!	----------------------------------------------
!	theta - Derivatives of the diffusion coefficients
    dKthetthet_dthet = 0.0

    ! TODO: DEBUG!!!
!    WRITE(*,*) "r, th, phi, E = ", r, theta, phi, E
!    !WRITE(*,*) "K_par = ", K_parallel
!    WRITE(*,*) "K_rr = ", K_rr
!    WRITE(*,*) "    Krr2 = ", K_rr2
!    WRITE(*,*) "K_phph = ", K_phiphi
!    WRITE(*,*) "K_rph = ", K_rphi
!    WRITE(*,*) "    Krph2 = ", K_rphi2
!    WRITE(*,*) "K_thth = ", K_thetatheta
!    WRITE(*,*) "dKrr_dr = ", dK_rrdr
!    WRITE(*,*) "    delta_r = ", delta_r
!    WRITE(*,*) "dKrph_dr = ", dK_rphi_dr
!    WRITE(*,*) "vdCoeff = ", drift_coeff / FS
!    WRITE(*,*) "vdr = ", vdr
!    WRITE(*,*) "vdth = ", vdtheta
!    WRITE(*,*) "vdph = ", vdphi
    !WRITE(*,*) "|B| = ", Bmag

!	----------------------------------------------
!	Rest of the numerical coefficients

	AA = 2.0/r*K_rr + dK_rrdr - V_sw - vdr + 1.0/r/sin(theta)*Krphidphi
	CC = K_rphi/r/r/sin(theta) + dK_rphi_dr/r/sin(theta) - vdphi/r/sin(theta)

	MM = cos(theta)/sin(theta)*K_thetatheta/r/r - vdtheta/r + dKthetthet_dthet/r/r
	NN = SQRT(2.0*K_thetatheta)/r*sqrt(delta_T)

!	----------------------------------------------
!	Write some output

	IF (printer.EQ.0) THEN

!	Plasma profiles.txt
		!WRITE(600,"(8(ES18.8))") r, E, P, theta, V_sw, DivV_sw, Bmag/0.06, larmor

!	Diffusion Coefficients
		!WRITE(700,"(7(ES18.8))") r, E, P, theta, lambda_parallel, lambda_perpr, lambda_perptheta

		printer = printer + 1

	END IF

	RETURN
	END

!-----------------------------------------------------
	SUBROUTINE LIS(E_end,P,j_N,E_0)
		IMPLICIT NONE

	REAL :: E_end,P,j_N
	REAL :: j_LIS
	REAL :: theta, E_0, a, b, c, d


		P = SQRT(E_end*(E_end + 2.0*E_0))
 
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
	SUBROUTINE JOVIAN_LIS(E_end, E_0, j_JUP)
		IMPLICIT NONE

	REAL :: E_end,P,j_JUP
	REAL :: j_LIS
	REAL :: E_0, a, b, c, d, e, j1, j2

		P = SQRT(E_end*(E_end + 2.0*E_0))

	a = 5000.0
	b = 1000000000.0
	c = 3.0
	d = 0.6
	e = 0.5

	j1 = a*(E_end*1000.0)**(-1.5)
	j2 = b*(E_end*1000.0)**(-6.0)  	  
	j_LIS = c*(d*j1*e*j2)/(d*j1+e*j2)

		j_JUP = j_LIS/P/P


	RETURN
	END

!-----------------------------------------------------
