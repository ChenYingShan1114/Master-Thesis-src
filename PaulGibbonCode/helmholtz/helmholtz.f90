!     -------------------------------------------------
!     
!     Helmholtz wave solver for inhomogeneous plasma profiles
!     
!     ------------------------------------------------
      
program helmholtz
  implicit none
  integer, parameter :: nmax=10000
  complex, dimension(0:nmax) :: eps
  complex, dimension(nmax) :: ez, bz, ex, ey, by, ez2
  complex, dimension(nmax) :: eza, bza, exa, eya, bxa, bya
  real, dimension(nmax) :: rho, x, phi, nu
  real, dimension(nmax)::  bzr,ezr,byr,eyr,exr,bxr,epond_p,epond_s
  complex :: yi
  real :: L,nu0,nusolid, wp, pi,Rs,Rp,thetar,kx, ky, thetaopt
  real :: a0, theta, rho0, vosc, vte, x0, xsolid, xl, dx, aginz
  real :: xd, ls, v, phase_s, phase_p, eps_sol, eps_step, phase_cs, phase_cp
  real :: psi_s,psi_p
  integer :: nprof, i, n, ithetaopt, ithe

  pi = asin(1.0)*2
  a0 = 1.
  yi = (0.,1.)
  theta=60.
  rho0=4.
  nusolid=0.01
  L=0.02

  write (*,*) 'Input max. density, collision frequency, scalelength, theta'
  !      write(*,*) " Give vosc:"
  !      read (*,*) rho0,nusolid,L,theta
  thetar=pi*theta/180.
  vosc=1.
  wp = sqrt(rho0)
  vte = 1./511.

  !     collision frequency at critical
  nu0 = nusolid/rho0
  !     nu0 = 0.02
  nprof=1
  !     position of vacuum
  x0 =8


  if (nprof.eq.1) then
     xsolid = x0+L*rho0
  else 
     xsolid = x0
  endif
  xl = xsolid + 2
  n = 1000

  dx=xl/n

  !     compute Ginzburg res. abs. optimum angle
  aginz = sqrt(0.5/(L)**(2./3.))
  if (aginz.lt.1) then
     thetaopt = 180/pi*asin(aginz)
  else 
     thetaopt=45.
  endif
  !     round off
  !     ithetaopt = thetaopt
  ithetaopt = 0
  ithe=0

  do i=1,n
     x(i)=dx*i
     if (nprof.eq.1) then
        !     linear profile
        if (x(i).lt.x0) then
           rho(i) = 0.
           nu(i) = 0.
        else if (x(i).lt.xsolid) then
           xd = x(i)-x0
           rho(i) = xd/L
           nu(i) = nu0*xd/L
        else 
           rho(i) = rho0
           nu(i) = nu0
        endif
     else if (nprof.eq.2) then
        !     exponential
        xd = x(i)-x0
        if (x(i).lt.x0) then

           rho(i) = rho0*exp(xd/L)
           nu(i) = nusolid*exp(xd/L)*(1+vosc**2/vte**2)**(-1.5)
        else
           rho(i) = rho0
           ! correct for vosc
           v = vosc*exp(-.5*xd/L)
           nu(i) = nusolid*(1+v**2/vte**2)**(-1.5)
        endif

     else
        ! step - warning: present routine cannot handle discontinuous eps(x)!!
        xd = x(i)-x0
        if (x(i).lt.x0) then
           rho(i)=0.
           nu(i) = 0.
        else
           rho(i) = rho0
           nu(i) = nu0
        endif

     endif
     eps(i) = 1.-rho(i)/(1+yi*nu(i))
  end do

  eps_sol = eps(n)
  eps(n+1) = eps(n)

  write (*,*) 'scale-length:',L
  write (*,*) 'resonance opt:',thetaopt
  open(30,file='abs_theta.dat')


  !     loop over angle
  !      do ithe = 0,179,2

  call solveRs(a0,dx,theta,eps(0:n+1),ez,n)
  call solveRp(a0,dx,theta,eps(0:n+1),bz,n)

  !     output fields at res. abs. optimum
  if (ithe.eq.ithetaopt) then
     bz(n+1) = 0.
     ez(n+1) = 0.
     do i=1,n
        ex(i) = -sin(thetar)*bz(i)/eps(i)
        ey(i)= -yi*2*(bz(i+1)-bz(i))/(eps(i)+eps(i+1))/dx
        by(i) = yi*(ez(i+1)-ez(i))/dx
     end do

     ! Analytical solutions

     kx=cos(thetar)
     ky=sin(thetar)
     ls=1./wp/sqrt(1-cos(thetar)**2/wp**2)
     !   ls=1./(sqrt(1.-kxokp^2))

     phase_cs = 2*atan(-1./kx/ls)
     phase_cp = 2*atan(-1./kx/ls/eps_sol)
     phase_s = atan(-kx*ls)
     phase_p = atan(-kx*ls*eps_sol)

     write (*,*) 'phi_s, phi_p: ',phase_cs,phase_cp,phase_s,phase_p
     do i=1,n
        xd = x(i)-x0
        ez2(i)= yi*2*(by(i+1)-by(i))/(eps(i)+eps(i+1))/dx
        if (xd.lt.0) then
           eza(i) = a0*(cexp(yi*kx*xd) + cexp(-yi*kx*xd+yi*phase_cs))
           bza(i) = a0*(cexp(yi*kx*xd) + cexp(-yi*kx*xd+yi*phase_cp))
           bzr(i) = -2*a0*sin(kx*xd+phase_p)
           ezr(i) = 2*a0*sin(kx*xd+phase_s)

       else
           eza(i) = a0*(1+cexp(yi*phase_cs))*exp(-xd/ls)
           bza(i) = a0*(1+cexp(yi*phase_cp))*exp(-xd/ls)
           bzr(i) = -2*a0*sin(phase_p)*exp(-xd/ls)
           ezr(i) = 2*a0*sin(phase_s)*exp(-xd/ls)

        endif
     end do

! derived fields

     bza(n+1) = bza(n)
     eza(n+1) = eza(n)

     do i=1,n
        xd = x(i)-x0
        if (xd.lt.0) then
           eps_step=1.
           eya(i)= -yi*(bza(i)-bza(i-1))/(eps_step)/dx
           psi_p = kx*xd+phase_p
           psi_s = kx*xd+phase_s
           eyr(i) = 2*a0*kx*cos(psi_p)
           exr(i) = 2*a0*ky*sin(psi_p)
           byr(i) = 2*a0*kx*cos(psi_s)
           bxr(i) = 2*a0*ky*sin(psi_s)
           epond_p(i) = -2*a0**2*kx*sin(2*psi_p)
           epond_s(i) = 2*a0**2*kx*sin(2*psi_s)
        else
           eps_step=1.-wp**2
           eya(i)= -yi*(bza(i+1)-bza(i))/(eps_step)/dx
           eyr(i) = -2*a0/eps_step/ls*sin(phase_p)*exp(-xd/ls)
           exr(i) = 2*a0*ky/eps_step*sin(phase_p)*exp(-xd/ls)
           byr(i) = -2*a0*sin(phase_s)/ls*exp(-xd/ls)
           bxr(i) = 2*a0*ky*sin(phase_s)*exp(-xd/ls)
           epond_p(i) = 4*a0**2/eps_step*sin(phase_p)**2/ls*exp(-2*xd/ls)
           epond_s(i) = -4*a0**2*sin(phase_s)**2/ls*exp(-2*xd/ls)
        endif

        exa(i) = -sin(thetar)*bza(i)/eps_step
        bya(i) = yi*(eza(i+1)-eza(i))/dx
        bxa(i) = sin(thetar)*ez(i)

     end do

     !  write out fields in units of c/wp
     open(20,file='field.dat')
     write(20,'(10(1pe15.5))') ((x(i)-x0)*wp,rho(i) &
          ,abs(ez(i)),abs(by(i)**2),abs(ex(i)),abs(ey(i)) &
          ,abs(ez(i))**2,abs(bz(i))**2,abs(ex(i))**2 &
          ,abs(ey(i))**2 &
          ,i=1,n)
     close(20)



     open(21,file='field_ana.dat')
     write(21,'(10(1pe15.5))') ((x(i)-x0)*wp,ezr(i)**2 &
          ,bzr(i)**2 &
          ,exr(i)**2,eyr(i)**2,bxr(i)**2,byr(i)**2 &
          ,epond_p(i),epond_s(i),real(eps(i)),i=1,n)
     close(21)

     !     reflection coeff

     !     for s-pol have E(1) = E_i + rE_r
     Rs = abs(ez(1)/a0-1.)**2
     !     for p-pol have B(1) = B_i - rB_r
     Rp = abs(1.-bz(1)/a0)**2
     write (*,*) 's-pol'
     write (*,*) 'Rs: ',Rs,'  As:',1-Rs
     write (*,*) 'p-pol'
     write (*,*) 'I/Ap: ',1.4e18*vosc**2,1-Rp
  endif

  !     reflection coeff

  !     for s-pol have E(1) = E_i + rE_r
  Rs = abs(ez(1)/a0-1.)**2
  !     for p-pol have B(1) = B_i - rB_r
  Rp = abs(1.-bz(1)/a0)**2

  write(30,'(3(f15.5))') theta,1-Rs,1-Rp
  !      end do


  !      write (*,*) 's-pol'
  !      write (*,*) 'Reflectivity: ',Rs,'  Absorption:',1-Rs
  !      Rwkb = exp(-32./15.*nu0*L*cos(pi/180*theta)**5)
  !      write (*,*) 'WKB result:',Rwkb


  close(30)
end program helmholtz

! -------------------------------------------------


!  solves Helmholtz equation for s-polarised light
!  - electric field Ez

subroutine solveRs(a0,dx,theta,eps,z,n)

  complex :: alpha(n),beta(n),gamma(n),y(n),z(n),eps(0:n+1)
  complex :: yi
  pi = asin(1.0)*2

  yi = (0.,1.)
  s2th=sin(pi/180*theta)**2
  cth = cos(pi/180*theta)

  do i=1,n
     !  coefficients
     y(i) = (0.,0.)
     alpha(i)=1
     beta(i)=-2 + dx**2*(eps(i)-s2th)
     gamma(i)=1
  end do

  !  BCs
  y(1) = 2*yi*a0*sin(dx*cth)
  beta(1) = beta(1) + cexp(yi*dx*cth)

  z(n)=(0.,0.)

  call trisolve(alpha,beta,gamma,y,z,n-1,n)

end subroutine solveRs

! -------------------------------------------------

!  solves Helmholtz equation for p-polarised light
!  - magnetic field Bz

subroutine solveRp(a0,dx,theta,eps,z,n)

  complex, dimension(n) :: alpha, beta, gamma, y, z
  complex :: eps(0:n+1)
  complex :: yi,deps
  pi = asin(1.0)*2
  yi = (0.,1.)
  s2th=sin(pi/180*theta)**2
  cth = cos(pi/180*theta)
  !      write (*,*) cth,s2th
  eps(0) = eps(1)
  eps(n+1) = eps(n)

  do i=1,n
     deps = (eps(i+1)-eps(i-1))/4.
     !  coefficients
     y(i) = (0.,0.)
     z(i) = (0.,0.)
     alpha(i) = 1. + deps/eps(i)
     beta(i) = -2. + dx**2*(eps(i)-s2th)
     gamma(i) = 1. - deps/eps(i)
  end do

  y(1) = 2*yi*a0*sin(dx*cth)
  beta(1) = beta(1) + cexp(yi*dx*cth)

  z(n)=(0.,0.)

  call trisolve(alpha,beta,gamma,y,z,n-1,n)


end subroutine solveRp

! -------------------------------------------------
!
!   Triadiagonal matrix solver.
!
!   Solves equation system:
!
!        alpha_i x_(i-1) + beta_i x_i + gamma_i x_(i+1) = y_i
!
! ------------------------------------------------

subroutine trisolve(alpha,beta,gamma,y,x,n,nmax)
  complex, dimension(nmax) :: alpha, beta, gamma,y,x, q

  q(1) = beta(1)
  x(1) = y(1)/q(1)

  !  forward elimination

  do i = 2,n
     q(i) = beta(i) - alpha(i)*gamma(i-1)/q(i-1)
     x(i) = (y(i) - alpha(i)*x(i-1))/q(i)
  end do

  !  back substitution

  do i=n-1,1,-1
     x(i) = x(i) - gamma(i)/q(i)*x(i+1)
  end do
end subroutine trisolve



!r1prnt
subroutine r1prnt(a,r1,m)
  character*6 a
  dimension r1(m)
  write(6,10) a,(r1(i),i=1,m)
10 format(1x/1x,a6/(1x,6(1pe12.3)))
  return
end subroutine r1prnt


!     c1prnt
subroutine c1prnt(a,c1,m)
  complex c1
  real a
  dimension c1(m)
  dimension c1r(5000),c1i(5000)
  do  i=1,m
     c1r(i)=real(c1(i))
     c1i(i)=aimag(c1(i))
  end do
  write(6,10) a,(c1r(i),i=1,m)
  write(6,11) a,(c1i(i),i=1,m)
10 format(1x/1x,a6,' (real part)'/(1x,6(1pe12.3)))
11 format(1x/1x,a6,' (imag part)'/(1x,6(1pe12.3)))
  return
end subroutine c1prnt
