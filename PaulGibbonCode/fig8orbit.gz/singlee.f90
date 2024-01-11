program singlee

parameter(n=100000)

!  Declare field and particle variables
real :: Ey, Bz, a0, px, py, x, y, q_over_m, gamma
real :: trun, delta_t, pi

!  intermediate variables
real :: beta, px_minus, py_minus, px_plus, px_dash, py_plus, t, s, gamma_n

pi = 2.*asin(1.) 
a0 = 1.0                   ! laser amplitude 
trun = 3*2*pi              ! run time (normalised laser periods)
delta_t = 0.1              ! timestep
nstep = trun/delta_t       ! # integration steps

!  Open a file for storing the particle orbit

open(30,file='orbit.dat')

x = 0.             ! initial position of electron
y = 0.
px = 0.            ! initial momentum of electron
py = 0.
q_over_m = -1.0    ! charge:mass ratio

time = 0.


! -----------------
! **  Main loop  **
! -----------------

do while (time.le.trun)

!  update laser field

    phase = time - x
    Ey = a0*sin(phase)
    Ex = 0.
    Bz = a0*sin(phase)

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

!  write out positions & momenta to file
    
    write (30,'(4f10.5)') x,y,px,py

!  increment time variable

    time = time + delta_t

end do


!  close files

close(30)

end
