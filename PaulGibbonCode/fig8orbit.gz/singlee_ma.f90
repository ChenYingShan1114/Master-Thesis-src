program singlee

parameter(n=100000)

!  Declare field and particle variables

real :: Ey, Bz, a0, px, py, x, y 

pi = 2.*asin(1.) 
a0 = 3.0                   ! laser amplitude 
tpulse =500
!tpulse = 0.
sigma = 20. 		  ! pulse width
trun = 500+2*pi              ! run time (normalised laser periods)
delta_t = 0.2              ! timestep
nstep = trun/delta_t       ! # integration steps
gamma0 = sqrt(1 + a0**2/2)  ! Sarachik & Schappert gamma

iseed = 7877
iseed2 = 90235
nsample = 200	 	! # samples for statistical run

q_over_m = -1.0    ! charge:mass ratio

!  Open a file for storing the particle orbit
open(30,file='orbit.dat')

!  File for storing final momenta
open(40,file='momenta.dat')

!  Loop over test particles
ipart = 1

do while (ipart <= nsample)


if (nsample.eq.0) then
  x = 0.             ! initial position of electron
  y = 0.1
else
  x = 2*pi*ran(iseed)
  y = 3*sigma*ran(iseed2)
  write(6,*) 'x0, y0:',ipart,x,y
endif

! initial momentum of electron

if (tpulse.eq.0) then
!  lab frame sawtooth solution
  py = a0
!  px = a0**2/2   

!  ave. rest frame
  px = a0**2/4/gamma0  
   
else
	!  start electron from rest
  py = 0.
  px = 0.
endif



time = 0.

! -----------------
! **  Main loop  **
! -----------------

do while (time.le.trun)

!  update laser field

    phase = time - x

    if (tpulse.gt.0) then
!  time envelope - sin2 pulse, duration = tpulse
      if (time.le.tpulse) then
        tenv = sin(pi*time/tpulse)**2
      else 
        tenv = 0.
      endif
    else 
!  constant amplitude
      tenv = 1.
    endif

    if (sigma.gt.0) then
!  pulse envelope - Gaussian, 1/e width = sigma
      earg = min(20.,y**2/2./sigma**2)
      space_env = exp(-earg)
    else 
!  constant amplitude
      space_env = 1.
    endif

    Alaser = a0*tenv*space_env*cos(phase)
    Ey = a0*tenv*space_env*sin(phase)
    Ex = 0.
    Bz = a0*tenv*space_env*sin(phase)

!  Buneman 'particle pusher' algorithm:

    beta = 0.5*q_over_m * delta_t     ! constant absorbing charge, mass and dt

    px_minus = px + beta*Ex           ! 1/2-acceleration
    py_minus = py + beta*Ey

    gamma_n = sqrt(1 + px_minus**2 + py_minus**2)

    t = beta*Bz/gamma_n               ! intermediate variables
    s = 2*t/(1 + t**2)

    px_dash = px_minus + py_minus*t    ! 
    py_plus = py_minus - px_dash*s     ! rotation about B
    px_plus = px_dash  + py_plus*t     !

    px = px_plus + beta*Ex            ! 1/2-acceleration
    py = py_plus + beta*Ey

    gamma = sqrt(1 + px**2 + py**2)    ! new relativistic factor

! compute velocities

    vx = px/gamma
    vy = py/gamma

! update particle position

    x = x + vx*delta_t
    y = y + vy*delta_t

 if (ipart.eq.193) then
!  write out positions & momenta to file

! x-coord normalised to lab-frame period
!    write (30,'(6f10.4)') time,2*x/pi/a0**2,y,px,py,a0*tenv*cos(phase)

!  self-similar in ave. rest frame
    q =a0/2./gamma0
!    write (30,'(6f10.4)') time,2*x/q**2,y/2./q,px,py,a0*tenv*cos(phase)

!  raw data
    write (30,'(6f10.4)') time,x,y,px,py,Alaser
endif

!  increment time variable

    time = time + delta_t

end do

!  write out final momenta/kinetic energy vs angle
 theta_code = 180./pi*atan(py/px)
! theta_theory = 180./pi*atan(sqrt(2./(gamma-1.)))
 write (40,'(2f10.5)') (gamma-1.0)*511.,theta_code
 ipart = ipart + 1
end do

! write out theoretical curve
open(50,file='theory.dat')
write(50,'(2f10.5)') ((g-1)*511.,180./pi*atan(sqrt(2./(g-1.))),g=1.01,2.,.01)

!  close files

close(30)
close(40)
close(50)
end
