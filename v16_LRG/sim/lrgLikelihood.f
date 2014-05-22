! Calculates the log-likelihood from SDSS LRG data as described in Appendix A.4 of
! Tegmark, Eisenstein, Strauss, Weinberg et al 2006, astro-ph/0608632
! COMPILATION: f77 compute_lrg_likelihood.f
! If you have questions about this code, please ask Max Tegmark, tegmark@mit.edu.

      subroutine compute_lrg_likelihood
        implicit none
        include 'lrg_likelihood_common.f'
        real*8 result, p(9), computeChisqf
        call initializef
        
        p(1) = 0        ! Okega_k
        p(2) = 0.762          ! Omega_Lambda
        p(3) = 0.1045         ! omega_c = h^2 Omega_c
        p(4) = 0.02233  ! omega_b  = h^2 Omega_b
        p(5) = 0.951          ! n_s
        p(6) = 0.6845         ! A_s
        p(7) = 0        ! alpha = dn/dlnk        
        p(8) = 1.908          ! galaxy bias
        p(9) = 30.81          ! nonlinear correction Q_nl
      
        result = computeChisqf(9, p)
        write(*,'(a10, f12.8)') 'chisq', result
      end


      subroutine initializef
        logical initialized
        save initialized
        data initialized/.false./
        if (.not.initialized) then
           print *,'Initializing SDSS likelihood stuff...'
           call initialize_sdss
           initialized = .true.
        end if
      end

      real*8 FUNCTION computechisqf(n ,array)
        integer i
        include 'lrg_likelihood_common.f'
        real*8 a, chi2, lnL
        real*8 array(*), kscaled(maxbands), power(maxbands)

!       do i=1,n
!          write(*,'(I5 F10.4)') I, array(i)
!       enddo

        ! Compute the theoretical P(k): 
        call compute_scaling_factor_from_p(array, a)
!       print *,'Scaling factor a =',a
        do i=1,kbands
          kscaled(i) = kvec(i)*a 
        end do    
        call computeP(array, kbands, kscaled, power)

        do i=1,kbands
          power(i) = power(i)/a**3
        end do

        ! Compute chi2:
        call compute_sdss_chi2(power,chi2)
        computechisqf = chi2
        RETURN
      END

      
      program demo
      implicit none
      include 'lrg_likelihood_common.f'
      integer i, idx
      real*8 kscaled(maxbands), power(maxbands), p(9), a, lnL, chi2
      character*20 args(9)
      character*50 outFile
      logical initialized
      save initialized
      data initialized/.false./
      if (.not.initialized) then
c         print *,'Initializing SDSS likelihood stuff...'
        call initialize_sdss
        initialized = .true.
      end if

      ! read in p
      do idx = 1,9
        call getarg(idx, args(idx))
        read (args(idx),*) p(idx)
      end do
      call getarg(10, outFile)

      ! First specify which model we wish to test.
      ! This example corresponds to the best-fit WMAP3 model from Table 5 of Spergel et al, astro-ph/0603449:
c       p(1) = 0      ! Okega_k
c       p(2) = 0.762        ! Omega_Lambda
c       p(3) = 0.1045       ! omega_c = h^2 Omega_c
c       p(4) = 0.02233      ! omega_b  = h^2 Omega_b
c       p(5) = 0.951        ! n_s
c       p(6) = 0.6845       ! A_s
c       p(7) = 0      ! alpha = dn/dlnk            
c       p(8) = 1.908        ! galaxy bias
c       p(9) = 30.81        ! nonlinear correction Q_nl
 
      ! Compute the theoretical P(k): 
      call compute_scaling_factor_from_p(p,a)
c       print *,'Scaling factor a =',a
      do i=1,kbands
        kscaled(i) = kvec(i)*a 
      end do            
      call computeP(p,kbands,kscaled,power)     ! Replace this by your own P(k) routine based on CMBfast, CAMB, CMBeasy, etc
                                    ! Simply neglecting the dewiggling isn't that bad an approximation.
      do i=1,kbands
        power(i) = power(i)/a**3
      end do            
      ! Note: reanalyzing the galaxy data with the flat Om=0.25 model replaced by the one 
      ! currently being tested would shift the measured P(k) curve up to the left if a>1,
      ! so we instead shift the theory curve down to the right.
      ! The measured curve would shift as k->k/a, P->P*a**3, so we 
      ! shift the theory curve as k->k*a,. P->P/a**3.

      ! Compute chi2:
      call compute_sdss_chi2(power,chi2)
      lnL = -0.5d0*chi2

c       write(*,'(a27,f12.8)') 'Chisq', chi2
c       write(*,'(a27,f12.8)') 'SDSS LRG log-likelihood is ',lnL
c       write(*,'(a27,f12.8)') 'The answer should be       ',-9.59273150d0

      ! write likelihood to file
      open(20, FILE=outFile, STATUS='REPLACE')
      write(20, *) lnL
      close(20)
      return
      end
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! ROUTINES DEALING WITH SDSS LRG DATA !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initialize_sdss
      ! Load kvec, obs, sdev, W into common block
      implicit none
      include 'lrg_likelihood_common.f'
      integer i, j, lnblnk
      real*8 sdev(maxbands), tmp(5)
      character*80 fname, s
      
      fname = 'sdss_lrg_kbands.txt'
      call LoadVector(fname,maxbands,kbands,kvec)   ! So we know the k-values at which to compute P(k) 
c       print *,kbands,' k-bands loaded from ',fname(1:lnblnk(fname))

      fname = 'sdss_lrg_measurements.txt'  ! Table 1 from astro-ph/0608632
      open(2,file=fname,status='old')
      read(2,*) s
      read(2,*) s
      bands = 0
555   read(2,*,end=666) (tmp(i),i=1,5)
      bands = bands + 1
      if (bands.gt.maxbands) stop 'DEATH ERROR: MAXBANDS TOO SMALL'
      obs(bands)  = tmp(4)
      sdev(bands) = tmp(5)
      goto 555
666   close(2)
c       print *,bands,' measurements loaded from ',fname(1:lnblnk(fname))
      bands = 20              ! Recommended. Corresponds to kmax ~ 0.2h/Mpc.
      ! set bands = 12        ! Use only small k-values if you're not applying a nonlinear correction
c       print *,bands,' measurements will be used.'

      fname = 'sdss_lrg_windows.txt'
      call LoadMatrix(maxbands,bands,kbands,W,fname)
c       print *,bands,' window fns loaded from ',fname(1:lnblnk(fname))

      ! Prewhiten so that WLOG N=I:
      do i=1,bands
        obs(i) = obs(i)/sdev(i)
        do j=1,kbands
          W(i,j) = W(i,j)/sdev(i)
        end do
      end do
      ! Now chi2 = (obs-W p)^t (obs-W p) = |obs-Wp|^2.

      return
      end
      
      subroutine compute_sdss_chi2(power,chi2)
      ! Data points are decorrelated and prewhitened with covariance natrix = I, 
      ! so chi2 = (obs-W p)^t (obs-W p) = |obs-Wp|^2.
      implicit none
      include 'lrg_likelihood_common.f'
      integer i
      real*8 power(maxbands), Wp(maxbands), chi2
      call vecmatmul(maxbands,bands,kbands,W,power,Wp)
      chi2 = 0.d0
      do i=1,bands
        chi2 = chi2 + (obs(i) - Wp(i))**2
      end do
      return
      end
      
      subroutine vecmatmul(mp,m,n,A,b,c)
      ! Matrix multiplication c=Ab of the m x n matrix A
      ! and the vector b.
      ! The physical (declared) number of rows are mp.
      implicit none
      integer mp,m,n,i,j      
      real*8 A(mp,n), b(n), c(n), sum
      do i=1,m
        sum = 0.
        do j=1,n
          sum = sum + A(i,j)*b(j)
        end do
        c(i) = sum
      end do
      return
      end
      
      subroutine LoadVector(f,nmax,n,x)   
      implicit none
      integer nmax,n,lnblnk
      real*8 x(nmax), r
      character*80 f
      n = 0
      open(2,file=f,status='old')
555   read(2,*,end=666) r
      n = n + 1
c       if (n.gt.nmax) pause 'n>nmax in LoadVector'
      x(n) = r
      goto 555
666   close(2)
      !print *,n,' numbers read from ',f(1:lnblnk(f))
      return
      end

      subroutine LoadMatrix(np,n,m,A,f)
      ! Reads the n x m matrix A from the file named f
      implicit none
      integer np,n,m,i,j,lnblnk
      real*8 A(np,m)
      character*80 f
      open(2,file=f,status='old')
      do i=1,n
        read(2,*) (A(i,j),j=1,m)
      end do      
      close(2)
      !print *,'Read matrix ',f(1:lnblnk(f))
      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! ROUTINES DEALING WITH GEOMETRIC SCALING FACTOR  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine compute_scaling_factor(Ok,Ol,w,a)
      ! a = dV for z=0.35 relative to its value for flat Om=0.25 model.
      ! This is the factor by which the P(k) measurement would shift sideways relative to what we got for this fiducial flat model.
      ! * a = (a_angular**2 * a_radial)**(1/3)
      ! * a_angular = comoving distance to z=0.35 in Mpc/h relative to its value for flat Om=0.25 model
      !     dA = (c/H)*eta = (2997.92458 Mpc/h)*eta, so we care only about eta scaling, not h scaling.
      !     For flat w=-1 models, a ~ (Om/0.25)**(-0.065)
      !     For the LRG mean redshift z=0.35, the power law fit dA(z,Om= 0.3253 (Om/0.25)^{-0.065}c H_0^{-1} is quite good within our range of interest,
        !     accurate to within about 0.1% for 0.2<Om<0.3. 
      ! * a_radial = 1/H(z) relative to its value for flat Om=0.25 model
      implicit none
      real*8 Or, Om, Ok, Ol, w, Ok0, Om0, Ol0, w0, z, eta, eta0, 
     +       Hrelinv, Hrelinv0, tmp
      real*8 a, a_radial, a_angular
      !Or   = 0.0000415996*(T_0/2.726)**4 / h**2
      Or    = 0   ! Radiation density totally negligible at  z < 0.35
      Om    = 1-Ok-Ol-Or
      z     = 0.35
      Hrelinv     = 1/sqrt(Ol*(1+z)**(3*(1+w)) + 
     +              Ok*(1+z)**2 + Om*(1+z)**3 + Or*(1+z)**4)
      call compute_z_eta(Or,Ok,Ol,w,z,eta)
      tmp = sqrt(abs(Ok))
      if (Ok.lt.-1.d-6) eta = sin(tmp*eta)/tmp
      if (Ok.gt.1d-6)   eta = (exp(tmp*eta)-exp(-tmp*eta))/(2*tmp) ! sinh(tmp*eta)/tmp    
      Ok0   = 0
      Ol0   = 0.75
      w0    = -1
      Om0   = 1-Ok0-Ol0-Or
      call compute_z_eta(Or,Ok0,Ol0,w0,z,eta0)
      Hrelinv0= 1/sqrt(Ol0*(1+z)**(3*(1+w0)) + Ok0*(1+z)**2 +
     +          Om0*(1+z)**3 + Or*(1+z)**4)
!     a_angular = (Om/0.25)**(-0.065) * (-w*Otot)**0.14 ! Approximation based on Taylor expansion
      a_angular = eta/eta0
      a_radial= Hrelinv/Hrelinv0
      a     =  (a_angular**2 * a_radial)**(1/3.d0)
      !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial,a_angular,(Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14
      !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial**(2/3.d0),a_angular**(4/3.d0),((Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14)**(4/3.d0)
      end
            
      subroutine compute_scaling_factor_from_p(p,a)
      !  Input parameters:
      !  1: Ok = Okega_k
      !  2: Ol = Omega_Lambda
      !  3: oc = omega_c = h^2 Omega_c
      !  4: ob = omega_b = h^2 Omega_b
      !  5: ns = n_s
      !  6: As = A_s
      !  7: al = alpha = dn/dlnk          
      implicit none
      real*8 p(7), a, Ok, Ol, w
      Ok = p(1)
      Ol = p(2)
      w  = -1
      call compute_scaling_factor(Ok,Ol,w,a)    
      return
      end

      subroutine eta_demo
      implicit none
      real*8 Or, Ok, Ol, w, h, z, eta
      h  = 0.7
      Ok = 0
      Ol = 0.7
      Or = 0.0000416/h**2 
      w  = -1
      z  = 1090
      call compute_z_eta(Or,Ok,Ol,w,z,eta)
      print *,'eta.............',eta
      print *,'dlss in Gpc.....',(2.99792458/h)*eta
      end
       
      logical function nobigbang2(Ok,Ol,w)
      ! Test if we're in the forbidden zone where the integrand blows up
      ! (where's H^2 < 0 for some a between 0 and 1).
      ! The function f(a) = Omega_m + Omega_k*a + Omega_l*a**(-3*w)
      ! can have at most one local minimum, at (Ok/(3*w*Ol))**(-1/(1+3*w)), 
      ! so simply check if f(a)<0 there or at the endpoints a=0, a=1.
      ! g(0) = Omega_m - Omega_l*a**(-3*w) < 0 if w > 0 & Omega_k > 1
      !                                     or if w = 0 & Omega_l < 1       
      ! g(1) = Omega_m + Omega_k + Omega_l = 1 > 0
      implicit none
      real*8 Ok, Ol, w, Om, tmp, a, epsilon
      integer failure
      failure = 0
      epsilon = 0
      !epsilon = 0.04  ! Numerical integration fails even before H^2 goes negative.
      Om = 1.d0 - Ok - Ol
      if (w*Ol.ne.0) then
        tmp = Ok/(3*w*Ol)
        if ((tmp.gt.0).and.(1+3*w.ne.0)) then ! f'(0)=0 for some a>0
          a = tmp**(-1/(1+3*w))
          if (a.lt.1) then
            if (Om + Ok*a + Ol*a**(-3*w).lt.epsilon) failure = 1
          end if
        end if
      end if
      if ((w.eq.0).and.(Ok.gt.1)) failure = 2
      if ((w.gt.0).and.(Ol.lt.0)) failure = 3
      nobigbang2 = (failure.gt.0)
      if (failure.gt.0) print *,'Big Bang failure mode ',failure
      return
      end

      real*8 function eta_integrand(a)
      implicit none
      real*8 Or, Ok, Ox, w
      common/eta/Or, Ok, Ox, w
      real*8 a, Om
      ! eta = int (H0/H)dz = int (H0/H)(1+z)dln(1+z) = int (H0/H)/a dlna = int (H0/H)/a^2 da = 
      ! Integrand = (H0/H)/a^2
      ! (H/H0)**2 = Ox*a**(-3*(1+w)) + Ok/a**2 + Om/a**3 + Or/a**4 
      if (a.eq.0.d0) then 
        eta_integrand = 0.d0
      else
        Om = 1.d0 - Or - Ok - Ox
        eta_integrand = 1.d0/sqrt(Ox*a**(1-3*w) + Ok*a**2 + Om*a + Or)
      end if
      return
      end
            
      subroutine eta_z_integral(Omega_r,Omega_k,Omega_x,w_eos,z,eta)
      ! Computes eta as a function
      ! of the curvature Omega_k, the dark energy density Omega_x
      ! and its equation of state w.
      implicit none
      real*8 Or, Ok, Ox, w
      common/eta/Or, Ok, Ox, w
      real*8 Omega_r, Omega_k,Omega_x,w_eos, z, eta, epsabs, epsrel, 
     +       amin, amax
      external eta_integrand
      Or = Omega_r
      Ok = Omega_k
      Ox = Omega_x
      w  = w_eos
      epsabs  = 0
        epsrel  = 1.d-10
      amin  = 1/(1+z)
      amax  = 1   
      call qromb2(eta_integrand,amin,amax,epsabs,epsrel,eta)
      return
      end
      
      subroutine compute_z_eta(Or,Ok,Ox,w,z,eta)
      ! Computes the conformal distance eta(z)
      implicit none
      real*8 Or, Ok, Ox, w, z, eta, Ht0
      logical nobigbang2
      if (nobigbang2(Ok,Ox,w)) then
        print *,'No big bang, so eta undefined if z>zmax.'
        eta = 99 
      else
        call eta_z_integral(Or,Ok,Ox,w,z,eta) 
        ! print *,'Or, Ok, Ox, w, z, H_0 t_0...',Or, Ok, Ox, w, eta
      end if
      return
      end
      
      SUBROUTINE qromb2(func,a,b,epsabs,epsrel,ss)
      ! The numerical recipes routine, but modified so that is decides
      ! it's done when either the relative OR the absolute accuracy has been attained.
      ! The old version used relative errors only, so it always failed when
      ! when the integrand was near zero.
      ! epsabs = epsrel = 1e-6 are canonical choices.
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,epsabs,epsrel
      EXTERNAL func
      PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.epsrel*abs(ss)) return
          if (abs(dss).le.epsabs) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      print *,'Too many steps in qromb'
      END

      SUBROUTINE polint(xa,ya,n,x,y,dy) ! From Numerical Recipes
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c           if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      SUBROUTINE trapzd(func,a,b,s,n) ! From Numerical Recipes
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! ROUTINES DEALING WITH THEORETICAL P(k) !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine computeP(p,kbands,kvec,power)
      !  Input parameters:
      !  1: Ok = Okega_k
      !  2: Ol = Omega_Lambda
      !  3: oc = omega_c = h^2 Omega_c
      !  4: ob = omega_b  = h^2 Omega_b
      !  5: ns = n_s
      !  6: As = A_s
      !  7: al = alpha = dn/dlnk          
      !  8: b  = galaxy bias
      !  9: Qnl= nonlinear correction Q_nl
      implicit none
      integer kbands, i
      real*8 p(9), kvec(kbands), power(kbands), b, Qnl, kstar
      b   = p(8)
      Qnl = p(9)
      call compute_kstar(p,kstar)
      call compute_EisensteinHu_dewiggled_power(p,kstar,kbands,kvec, 
     +                                          power)
      !print *,'Applying bias b=',b,' and Qnl-correction with Qnl=',Qnl
      do i=1,kbands
        power(i) = power(i) * b**2 * (1+Qnl*kvec(i)**2)/(1+1.4*kvec(i))
      end do
      return
      end

      subroutine compute_kstar(p,kstar)
      ! Computes the dewiggling parameter kstar using equations (12) & (13) in Eisenstein, Seo & White astro-ph 0604361
      implicit none
      real*8 p(7), kstar, As, Ok, Ol, Om, Okz, Olz, Omz, z, E, g, f,
     +       D, s
      Ok    = p(1)
      Ol    = p(2)
      !oc   = p(3)
      !ob   = p(4)
      !ns   = p(5)
      As    = p(6)
      !al   = p(7)
      Om    = 1-Ol-Ok
      z     = 0.35      ! Effective redshift of LRG's
      E     = Ol + Ok*(1+z)**2 + Om*(1+z)**3
      Olz   = Ol/E
      Okz   = Ok*(1+z)**2/E
      Omz   = Om*(1+z)**3/E
      g     = 2.5d0*Omz/(Omz**(4/7.d0)-Olz+(1+Omz/2)*(1+Olz/70)) ! g = D/a, The Carroll et al (1992) approximation
      f     = -Omz/2 - 1 + Olz + 2.5d0*Omz/g               ! f = dlnD/dlna; formula from Carroll, Press & Turner via Andrew Hamilton
      D     = g/(1+z)
      s     = 12.5 * sqrt(As/0.6841) * D * (1+f)**(1/3.d0)  ! Mpc/h
      kstar = 1/s
      !print *,'kstar=',kstar
      return
      end

      subroutine compute_EisensteinHu_dewiggled_power(p,kstar,bands,
     +                                                kvec,power)
      ! Uses the Eisenstein & Hu fitting formulas to compute P(k) rapidly.
      implicit none
      integer bands, i
      real*8 p(7), kstar, kvec(bands), power(bands)
      do i=1,bands
        call EH_computeP_dewiggled(p,kstar,kvec(i),power(i))
      end do
      end

      subroutine EH_computeP_dewiggled(p,kstar,k,power)
      ! Uses the Eisenstein & Hu fitting formulas to compute P(k) rapidly.
      ! k is in h/Mpc, *not* 1/Mpc.
      ! I can accelerate this routine by doing most of the calcs (incl norm call) only 
      ! ininitialization, storing results in common block. But it's so fast anyway
      ! that I don't really care.
      !  Input parameters:
      !  1: Ok = Okega_k
      !  2: Ol = Omega_Lambda
      !  3: oc = omega_c = h^2 Omega_c
      !  4: ob = omega_b  = h^2 Omega_b
      !  5: ns = n_s
      !  6: As = A_s
      !  7: al = alpha = dn/dlnk          
      !  kstar = scale of Gaussian suppression of baryon wiggles, also given in h/Mpc.
      !  kstar >= 1e5 gives identical result to EH_computeP, kstar <= 1e-5 gives no wiggles at all.
!     implicit real*8 (a-h,k,o-z)
      implicit none
      integer i
      real*8 k, kstar, kk, p(7), lastp(7), h, T_cmb, omhh, f_baryon
      real*8  Ok, Ol, Om, oc, ob, ns, As, al, Tnorm, cmbfastnorm
      real*8 transfer, tf_full, tf_baryon, tf_cdm, power, T0, T, pi
      real*8 transfer_wiggles, transfer_nowiggles, TF_nowiggles,
     +       wiggleweight, x
      parameter(T0=2.725d6) ! CMB temperature in uK
      logical initialized
      data lastp/666, 666, 666, 666, 666, 666, 666/
      save lastp
      pi = 4.d0*atan(1.d0)
      ! Input Cosmological Parameters
      Ok    = p(1)
      Ol    = p(2)
      oc    = p(3)
      ob    = p(4)
      ns    = p(5)
      As    = p(6)
      al    = p(7)
      Om    = 1.-Ol-Ok
      omhh  = oc + ob
      f_baryon= ob/omhh
      h     = sqrt(omhh/(1.-Ok-Ol))
      T_cmb = T0/1d6    ! CMB temperature in K
      kk = h*k          ! kk has units of of 1/Mpc
      initialized = .true.
      do i=1,4 ! Need to recompute transfer function if any of 1st four parameters changed
        if (p(i).ne.lastp(i)) initialized = .false.
      end do 
      if (.not.initialized) then
        do i=1,7
          lastp(i) = p(i)
        end do 
          call TFset_parameters(omhh,f_baryon,T_cmb)  
      endif       
      ! Eisenstein & Hu work in 1/Mpc, not h/Mpc
      call TFtransfer_function(kk,omhh,f_baryon,tf_full,tf_baryon,
     +                         tf_cdm)
      transfer_wiggles   = tf_full
      transfer_nowiggles = TF_nowiggles(kk,omhh,f_baryon,T_cmb)   
      if (kstar.le.1.d-5) then
        wiggleweight = 0
      else 
        x = k/kstar
        if (x.gt.10) then 
          wiggleweight = 0
        else
          wiggleweight = exp(-x**2/2)
        end if
      endif
      if (kstar.ge.1.d5) wiggleweight = 1
      transfer = wiggleweight*transfer_wiggles + 
     +           (1-wiggleweight)*transfer_nowiggles
      call ComputeT0(h,Om,Ol,Tnorm)
      T = Tnorm*transfer
      cmbfastnorm = 800*(pi*5./3.)**2/T0**2 ! Converts to our old normalization As=1
      ! The next three lines are the same as in cmbfast_compute_P.f: (writing kk in place of k)
      power = 2*pi**2 * (kk/0.05)**(ns - 1 + al*log(kk/0.05)/2) *
     +        kk * T**2 ! Now power has units of Mpc**3
      power = h**3 * power                       ! Now power has units of (Mpc/h)**3, normalized to A=1
      power = cmbfastnorm * As * power            
      return
      end

      subroutine ComputeT0(h,Om,Ol,Tnorm)
      ! Computes Tnorm, the CMBfast transfer function T(k) as k=0.
      implicit none
      real*8 h, Om, Ol, Tnorm, Tstar, growth
      Tstar = 3594770 ! This is empirical cmbfast T(0) for flat model with Om=1, h=1      
      growth = 2.5*Om/(Om**(4./7.)-Ol+(1+Om/2)*(1+Ol/70)) ! The Carroll et all approximation
      Tnorm = Tstar*growth/(Om*h**2)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!! HERE COMES THE EISENSTEIN & HU PACKAGE !!!!!!!!!!!!!!!!!!!!!!!!!
c     Driver and Transfer Function / Power Spectrum subroutines 
c            for Eisenstein & Hu astro-ph/9709112 
c       (This version does wiggles but not neutrinos.
c        the other does neutrinos but not wiggles.)

c   The following routines implement all of the fitting formulae in 
c   Eisenstein \& Hu (1997) 
c
c   Program TF_fit: sample driver
c   Subroutine TFset_parameters(): sets all the scalar parameters
c   Subroutine TFtransfer_function(): calculates various transfer functions
c   Functions TF_zerobaryon, TF_nowiggles, sound_horizon_fit, kpeak,
c       alpha_gamma: implement various scaling approximations of \S 4.2
c
c
c ------------------------ DRIVER ROUTINE --------------------------- 
c The following is a driver routines you might use. 
c Basically, the driver routine asks for Omega_0, the baryon fraction,
c the hubble constant, and the CMB temperature, calls TFset_parameters() to
c set all the parameters of the fit.  A loop over wavenumbers kmin to an
c inputed kmax sampled at numk per decade calls TFtransfer_function. 
c
c IMPORTANT: TFtransfer_function asks for wavenumbers in Mpc^{-1} so
c          multiply by hubble to convert from h Mpc^{-1}
c
c The latter returns values of the various transfer functions at the given 
c wavenumber which are output to the file "trans.dat" 
c
c Also included is an example of how to call the functions
c
c     TF_nowiggles:         shape approximation
c     TF_zerobaryon:      zero baryon TF 
c     k_peak:               approximate first peak location
c     sound_horizon_fit:  approximate sound horizon
c     alpha_gamma:          small scale modification to Gamma
c
c
c
c INPUT:  omega0 -- the matter density (baryons+CDM) in units of critical 
c       f_baryon -- the ratio of baryon density to matter density 
c       hubble -- the Hubble constant, in units of 100 km/s/Mpc
c       Tcmb -- the CMB temperature in Kelvin.  2.728(4) is COBE and is the default
c               reached by setting Tcmb=0.
c       kmax -- maximum k in  h Mpc^{-1}
c       numk -- number of k per decade 
c
c OUTPUT: The file "trans.dat" with columns
c
c      (1)  k -- wavenumber in h Mpc^{-1}
c      (2)  tf_full -- The full fitting formula, eq. (16), for the matter
c                 transfer function. 
c      (3)  tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
c      (4)  tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
c      (5)  tf_nowiggles -- An approximate form, eqs. (30)-(31), that fits
c                 only the non-oscillatory part of the transfer 
c                 function.  Appropriate only for low baryon fractions.
c      (6)  tf_zerobaryon -- The transfer function of the zero-baryon case,
c                 eq. (29); i.e. what would have occured were the
c                 baryons CDM instead.                
c
c     and the approximate values of k_peak,sound_horizon,alpha_gamma 
c     to stdout, a more accurate form of sound_horizon lives in
c
c     GLOBALVARIABLES: Various intermediate fit parameters are stored in 
c                        this common block for easy access.


c ----------------------------- DRIVER ------------------------------- 


       subroutine TFfit


        real*8    omega0,f_baryon,hubble,Tcmb,kmax,kmin, omhh
      real*8      k,tf_full,tf_baryon,tf_cdm,tf_nowiggles,tf_zerobaryon
      real*8  k_peak, sound_horizon_fit, alpha_gamma
      
      integer numk, i

c  cosmological parameters

      write(6,*) 'Omega_0,f_baryon,h,T_cmb?'
      read*,     omega0,f_baryon,hubble,Tcmb
      omhh = omega0*hubble*hubble

c  call routine to set fitting parameters

        call TFset_parameters(omhh, f_baryon, Tcmb)


c  loop over k, call subroutine and functions to calc TFs
 
      open(10,file='trans.dat')

      write(6,*) 'k_max (h Mpc^{-1}),#pts per decade (10,50)'
      read*,kmax,numk

      if (kmax.le.0) kmax=10.
      if (numk.le.0) numk=50

      kmin = 0.0001
      numk = numk*log10(kmax/kmin)

      do i=1,numk

       k=10.**(i*(log10(kmax/kmin)/numk))*kmin
         call TFtransfer_function(k*hubble,omhh,f_baryon,
     &     tf_full,tf_baryon,tf_cdm)
       write(10,50) k,tf_full,tf_baryon,tf_cdm,
     &               TF_nowiggles(k*hubble,omhh,f_baryon,Tcmb),
     &               TF_zerobaryon(k*hubble/omhh*(Tcmb/2.7)**2)


      end do

c  example of how to use the scaling functions

      write(6,*) 'Some useful approximate scalings:'

      write(6,10) k_peak(omhh,f_baryon)/hubble
10    FORMAT(1X,' First peak location (h Mpc^{-1}):  ',E13.4)
      write(6,20) sound_horizon_fit(omhh,f_baryon)*hubble
20    FORMAT(1X,' Approx. sound horizon (h^{-1} Mpc):',E13.4)
      write(6,30) alpha_gamma(omhh,f_baryon)
30    FORMAT(1X,' alpha_gamma:                       ',E13.4)


50    FORMAT(1X,7E13.5)

      end

c
c
c PART I:------------------- FITTING FORMULAE ROUTINES ----------------- 
c
c There are two routines and a set of functions.  
c   TFset_parameters() sets all the scalar parameters, while 
c   TFtransfer_function() calculates various transfer functions 
c
c Global variables -- We've left many of the intermediate results as
c global variables in case you wish to access them, e.g. by declaring
c them as a common block in your main program. 
c
c Note that all internal scales are in Mpc, without any Hubble constants! 
c

      subroutine TFset_parameters(omhh,f_baryon,Tcmb)

      real*8 y,omhh,obhh,Tcmb
      real*8 theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &           f_baryon,beta_node
      common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node

c Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
c Input omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c       f_baryon -- The fraction of baryons to CDM */
c       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
c           the default reached by inputing Tcmb=0 -- reset on output. */
c Output nothing, but set many global variables in common block 
c       GLOBALVARIABLES. You can access them yourself, if you want:
c
c     theta_cmb,  /* Tcmb in units of 2.7 K */ 
c     z_equality, /* Redshift of matter-radiation equality, real*8ly 1+z */
c     k_equality, /* Scale of equality, in Mpc^-1 */
c     z_drag,           /* Redshift of drag epoch */
c     R_drag,           /* Photon-baryon ratio at drag epoch */
c     R_equality, /* Photon-baryon ratio at equality epoch */
c     sound_horizon,    /* Sound horizon at drag epoch, in Mpc */
c     k_silk,           /* Silk damping scale, in Mpc^-1 */
c     alpha_c,    /* CDM suppression */
c     beta_c,           /* CDM log shift */
c     alpha_b,    /* Baryon suppression */
c     beta_b,           /* Baryon envelope shift */



c Are inputs reasonable?

      if (f_baryon.le.0) f_baryon=1.d-5
      if (Tcmb.le.0) Tcmb=2.728
        if (omhh.le.0.0) then
         write(6,*) 'TFset_parameters(): Illegal input'  
c          pause
      end if

!        if (hubble.gt.10.0) then
!        write(6,*) 'TFset_parameters(): WARNING, Hubble constant in 
!     &                    100km/s/Mpc desired'
!     end if

c Auxiliary variables
        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

c Main variables
        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.) 

          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1d0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291d0 * omhh**(0.251)/
     &           (1d0 + 0.659*omhh**(0.828)) * z_drag
 
        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000d0/(1d0 + z_drag) 
        R_equality = 31.5*obhh*theta_cmb**(-4.) 
     &           *1000d0/(1d0 + z_equality) 

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
     &          log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
     &       /(1.+sqrt(R_equality)))

        k_silk = 1.6*obhh**(0.52)*omhh**(0.73)* 
     &           (1d0 + (10.4*omhh)**(-0.95))

          alpha_c = ((46.9*omhh)**(0.670)*(1d0+(32.1*omhh)**(-0.532)))
          alpha_c = alpha_c**(-f_baryon) 
      alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1d0 + 
     &             (45.0*omhh)**(-0.582)))**(-f_baryon**3.)

    
          beta_c = 0.944/(1+(458.*omhh)**(-0.708))
          beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266)) 
     &      - 1d0)
        beta_c = 1./beta_c

          y = (1d0+z_equality)/(1d0+z_drag)
          alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)
     &          /(sqrt(1.+y)-1.)))
        alpha_b = 2.07*k_equality*sound_horizon*
     &            (1.+R_drag)**(-0.75)*alpha_b


        beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*
     &           sqrt((17.2*omhh)**2.+1d0)

        beta_node = 8.41*omhh**(0.435)

        return

        end


        subroutine TFtransfer_function(k,omhh,f_baryon,tf_full,
     &                tf_baryon,tf_cdm)

c  Calculate transfer function from the fitting parameters stored in
c  GLOBALVARIABLES.
c
c  Input: 
c      k -- wavenumber in Mpc^{-1}  
c        omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c        f_baryon -- The fraction of baryons to CDM */
c     
c  Output:
c      tf_full -- The full fitting formula, eq. (16), for the matter
c                 transfer function. 
c      tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
c      tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
c


      real*8 k,tf_full,tf_baryon,tf_cdm,q,ks, omhh, s_tilde
      real*8 theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &           f_baryon,beta_node, TF_pressureless
      common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node


c  Reasonable k?

        if (k.le.0) then
           write(6,*) 'TFtransfer_function(): Illegal k'
c            pause
        end if 


c  Auxiliary Variables

          q = k/13.41/k_equality
          ks = k*sound_horizon

c  Main Variables

              tf_cdm = 1./(1.+(ks/5.4)**4.)
          tf_cdm = tf_cdm*TF_pressureless(q,1.d0,beta_c) + 
     &           (1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)


            s_tilde = sound_horizon/(1.+(beta_node/ks)**3)**(1./3.) 
            tf_baryon = TF_pressureless(q,1.d0,1.d0)/
     +                                 (1.d0+(ks/5.2d0)**2)
            tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3)
     &                       *exp(-(k/k_silk)**(1.4))
            tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
          tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

         return

       end

c       auxiliary function: Pressureless TF

      real*8 function TF_pressureless(q,a,b)

        real*8 q,a,b
      
        TF_pressureless = Log(exp(1.)+1.8*b*q)
        TF_pressureless = TF_pressureless/(TF_pressureless + 
     &                      (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

        return

      end   
c
c
c
c
c
c
c PART II:------------------- Scaling Functions ROUTINES ----------------- 
c
c       omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c       f_baryon -- The fraction of baryons to CDM */
c
c
c     TF_zerobaryon:     
c       Input:  q = k/omhh * (Tcmb/2.7)**2    (k in Mpc^{-1})
c       Output: zero baryon TF Eq(29)
c     TF_nowiggles:      
c       Input:  k = wavenumber in Mpc^{-1}, omhh, f_baryon, Tcmb
c       Output: shape approximation TF  Eq(30-31)
c       Calls: TF_zerobaryon,sound_horizon_fit,alpha_gamma
c     sound_horizon_fit: 
c         Input:  omhh,f_baryon     
c       Output: approximate sound horizon in Mpc      
c     kpeak:               
c       Input:  omhh,f_baryon
c         Output: first peak location in Mpc
c       Calls:  sound_horizon_fit
c     alpha_gamma:         
c       Input: omhh,f_baryon
c       Output: effective small scale suppression


      real*8 function TF_zerobaryon(q)

        real*8 q
        TF_zerobaryon = log(2.0*exp(1.)+1.8*q)
        TF_zerobaryon = TF_zerobaryon/(TF_zerobaryon
     &                    +(14.2 + 731.0/(1+62.5*q))*q**2)

        return

      end

      real*8 function TF_nowiggles(k,omhh,f_baryon,Tcmb)

        real*8 k,omhh,f_baryon,q_eff,a,Tcmb
        real*8 sound_horizon_fit,alpha_gamma,TF_zerobaryon
      
          if (Tcmb.le.0) Tcmb=2.728
          a = alpha_gamma(omhh,f_baryon)
          q_eff = k/omhh*(Tcmb/2.7)**2
          q_eff = q_eff/(a+(1.-a)/
     &              (1.+(0.43*k*sound_horizon_fit(omhh,f_baryon))**4))

          TF_nowiggles = TF_zerobaryon(q_eff)
        
        return
        end


      real*8 function sound_horizon_fit(omhh,f_baryon)

       real*8 omhh,obhh,f_baryon
       obhh = f_baryon*omhh
         sound_horizon_fit = 44.5*log(9.83/omhh)
     &                      /sqrt(1.+10.0*obhh**(0.75))

      return
      end


      real*8 function k_peak(omhh,f_baryon)

       real*8 omhh,obhh,f_baryon,sound_horizon_fit
       obhh = f_baryon*omhh
         k_peak = 5.*3.14159/2.*(1.+0.217*omhh)/
     &             sound_horizon_fit(omhh,f_baryon)

       return
      end


      real*8 function alpha_gamma(omhh,f_baryon)

       real*8 omhh,f_baryon

         alpha_gamma = 1.-0.328*log(431.0*omhh)*f_baryon 
     &                + 0.38*log(22.3*omhh)*(f_baryon)**2
    
       return
      end 
