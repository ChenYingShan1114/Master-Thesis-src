      parameter(nrmax=1001,nzmax=20000,ntmax=10000,ncmax=10
     :         ,nzm=nzmax,nfmax=8192, maxspc=2
     :         ,nvxm=100
     :         ,pi=3.14159, c=2.998e8 )

      real me,mi,miome,n0
      complex a,yi

      parameter( yi=(0.,1.) )
      character ctime*12, gas*2
  
      real*8 rhoi, rhoatom, rhoe

c em and plasma arrays
      common /rgrid/
     :  a(0:nrmax,0:nzmax)
     :,rhoi(0:nrmax,0:nzmax,0:maxspc),rhoe(0:nrmax,0:nzmax)
     :,rhoatom(0:nrmax,0:nzmax)
     :,r(0:nrmax), xnu(0:nrmax), zdens(0:nzmax)
     :,phi(0:nrmax,0:nzmax), epond(0:nrmax,0:nzmax)
     :,Er(0:nrmax,0:nzmax)
     :,frion(0:nrmax)
     :,gamma(0:nrmax,0:nzm)

c  history diagnostics
      common /hist/
     : edena(0:ntmax), relem(0:ntmax), Ues(0:ntmax), wbk(0:ntmax)
     :,action(0:ntmax), rbeam(0:ntmax), H(0:ntmax), a2r0(0:ntmax)
     :,Pcon(ncmax), rcon(ncmax,0:ntmax), rhoer0(0:ntmax),rhoar0(0:ntmax)
     :,H2(0:ntmax),rho0r0(0:ntmax),rho1r0(0:ntmax),rho2r0(0:ntmax)
     :,a2ne(0:ntmax), aphi(0:ntmax)
     :,utot(0:ntmax), uth(0:ntmax),  udr(0:ntmax)
     :,esoth(0:ntmax), erh(0:ntmax), dene0(0:ntmax), denem(0:ntmax)


      common /physcon/
c  dimensionLESS
     : wp0, sigma, tau, t0, rl, zl, Zrmass, Zion, power, popc
     :,zedge, zend, zplas, timec, zbord, zfoff
     :,Ea0, wa0, rho0, rhoe0, zprobe, tprobe, omemin, omemax, zplen
c  laser
     :,a0, Rfoc, Zray, Dlens, fno, sigma0, xI0, ptw
     :,xlam, w0, trise, tdel, tp, fpol, trun
c  dimensionAL
     :,rfac, tfac, sfac, rplmax, den0, rthom
c  switches
     :,ilas, mabs, imovie, igrid, icycem, itav, izsurf
     :,melec, iout, mrel, mfoc, mioniz, mplas, nion
c  flags
     :,wbcm
c  diags
     :,ncon, isurf, ksurf, kprobe, ipskip, ipcon
     :,gas

      common/syst/
     :     wp,    n0,   xm1,   xm2,  xsol
     :,   abo,   asm,  ampl, xsol2,   tav,  gnon,  xgin,  gam0
     :,  erhb,    wl,    wr,  uth0
     :,xdebye,xdodr ,   xla,   rwp

      common /gricon/
     : nr, nz, nzp, nt, itime
     :,dt, dr, dz, dzp, rdr, dkx, dts, dtioniz, dtem
     :,igr, ihist, inorm
     :,ctime, nitem, nabs, xnumin, xnumax
     :,nro2,nfo2,dto2,isubm

      namelist /idynf/ 
     : wp0, sigma, tau, t0, rl, zl, Zrmass, dt, a0, Zion
     :,nr, nz, nzp, nt, igr, ihist, ilas, nitem, iout
     :, mplas, zplas, trun,ptw, rthom
     :, inprof, ipcon, izsurf, gas, isubm
     :,xlam, Rfoc, Dlens, fno, xi0, zplen, zbord, zfoff
     :,inorm, melec, imovie, igrid, isurf, ksurf,mrel, mfoc
     :,nabs, xnumin, xnumax, mabs, rplmax, fpol, den0
     :,pcon, ncon, mioniz, icycem, tprobe,omemin,omemax
     :,trise,xload,xm1,xm2,tp,tdel,nf,xsol,xsol2,Te




