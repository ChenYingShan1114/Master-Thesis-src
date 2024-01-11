!!
!! Program to solve Brunel's absorption model exactly
!! including relativistic and pump-depletion effects

program brunel

  implicit none
  integer, parameter :: nmax=200
  integer :: i
  real :: lambda, eta, a0, theta, pi, Ilam2, Intensity, ifix, fudge, eta_r
  real :: eta_theta(4), theta_deg, gm1
  pi = 2.*asin(1.0)
  write (*,*) 'Give theta (deg)'
  read (*,*)  theta
  lambda=0.4
  fudge=1.6
  theta = pi*theta/180.

! Intensity dependence
  open(30,file='abs_brunel_i.dat')

  do a0 = 0.02,10.,.1
     Ilam2 = 1.37e18*a0**2
     Intensity=Ilam2/lambda**2
     call solve_abs(theta,a0,fudge,eta)  ! Brunel model
     call solve_relabs(theta,a0,fudge,eta_r,gm1)  ! Fully relativistic VH
     write(30,'(f15.4,5(1pe12.3))') a0,Intensity,Ilam2,eta, eta_r,gm1*511.
  end do

  close(30)

! angular dependence
  write (*,*) " Angular dep"
  open(30,file='abs_brunel_theta.dat')

  do theta_deg = 0.,80.,1.
     theta=pi*theta_deg/180.
 !    a0=0.027
 !    call solve_relabs(theta,a0,fudge,eta_theta(1))  ! Fully relativistic VH
 !    a0=0.27
 !    call solve_relabs(theta,a0,fudge,eta_theta(2))  ! Fully relativistic VH
     a0=6.
     call solve_abs(theta,a0,fudge,eta_theta(3))  ! Fully relativistic VH

     write(30,'(f15.4,3(1pe12.3))') theta_deg,eta_theta(1),eta_theta(2),eta_theta(3)
  end do

  close(30)
end program brunel



subroutine solve_abs(theta,a0,fudge,eta)

implicit none
real, intent(out) :: eta
real, intent(in) :: theta, a0, fudge
real :: a,b,f,pi,alpha,eps,epsmax,rhs,beta
real :: eps_low, eps_high, eps_conv,trial
integer :: i

epsmax = 0.001  ! 1% error tolerance
eps=1.
trial = 0.1
! constants
pi = 2.*asin(1.0)
a = fudge/pi/a0*sin(theta)/cos(theta)
b = (a0*sin(theta))**2
!write (*,*) 'a0sin(theta) =',a0*sin(theta)
i=1
f = 1.+sqrt(1.0-trial)

do while (eps>epsmax .and. i<20)
   f = 2 - a*( sqrt(1. + b*f**2) - 1.)
   eta = 1.-(f-1)**2
   eps = abs(eta-trial)/eta
!   write (*,'(a2,i5,4(a10,f12.4))') 'i',i,'f ',f,'eta(n) ',trial,' eta(n+1)',eta,' error',eps
   trial = eta
!   f = 1.+sqrt(1.0-trial)
    i=i+1  
end do
if (i==20) then 
   write(*,*) "Not converged"
   write (*,'(5(a10,f12.4))') 'a0',a0,'f ',f,'eta(n) ',trial,' eta(n+1)',eta,' error',eps
else 
   write(*,*) "Done: eta=",eta," error",eps
endif
! limiting values
!beta=a0*(sin(theta))**3/cos(theta)/2./pi
f = (sqrt(1.+8*beta)-1.)/2/beta
!eps_low = 1-(f-1)**2 
eps_low = (4*a0/pi)*(sin(theta)**3/cos(theta))
alpha = (sin(theta))**2/cos(theta)
eps_high = 4*pi*alpha/(pi+alpha)**2
eps_conv = a*f*( sqrt(1. + b*f**2) - 1.)
!write (*,'(2(a4,f12.3),a10,1pe12.3,4(a15,1pe12.4))') 'theta',theta,'a0',a0,'Intensity',1.38e18*a0**2/0.4**2,'eta: ',eta,'Low limit: ',eps_low,' High limit: ',eps_high, 'conv',eps_conv
end subroutine solve_abs


subroutine solve_relabs(theta,a0,fudge, eta, gm1)

implicit none
real, intent(out) :: eta  ! absorption fraction
real, intent(out) :: gm1  ! Thot
real, intent(in) :: theta ! angle of incidence
real, intent(in) ::  a0  ! pump amplitude
real, intent(in) ::  fudge  ! empirical Brunel factor

real*8 :: a,b,c,f,pi,alpha,eps,epsmax,rhs,beta,gamma
real*8 :: eps_low, eps_high, eps_conv,trial
integer :: i

epsmax = 0.001  ! 1% error tolerance
eps=1.
trial = 0.1

! constants
pi = 2.*asin(1.0)
a = fudge/pi/a0*sin(theta)/cos(theta)
b = ( a0*sin(theta)*cos(theta) )**2
c = cos(theta)**2

!write (*,*) 'a0sin(theta) =',a0*sin(theta)
i=1
f = 1.+sqrt(1.0-trial)

do while (eps>epsmax .and. i<10)
   gamma = 1.+(sqrt(1.+b*f**2) - 1.)/c
   f = min(1.999,2. - a*(gamma-1.)/gamma)
   eta = 1.-(f-1.)**2
   eps = abs(eta-trial)/eta
!   write (*,'(6(a10,f12.6))') 'theta',theta,'a0',a0,'f ',f,'eta(n) ',trial,' eta(n+1)',eta,' error',eps
   trial = eta
!   f = 1.+sqrt(1.0-trial)
     i=i+1 
end do

if (i==10) then 
   write(*,*) "Not converged"
   write (*,'(5(a10,f12.4))') 'a0',a0,'f ',f,'eta(n) ',trial,' eta(n+1)',eta,' error',eps
else 
   write(*,*) "Done: eta=",eta," error",eps
endif
! Final K.E.
   gm1 = (sqrt(1.+b*f**2) - 1.)/c

end subroutine solve_relabs













