!
! Simulates the orbital evolution of a population
! of comets, originating in a Kuiper belt, using
! a simple 1D Markov-chain model.
!
! Cometary orbits are perturbed by a chain of planets.
! Comets are initially located just beyond outermost planet.
!
! Planetary masses are
!   M = M_small                         a < a_tran
!     = power-law function of a         a > a_tran
! where a is distance from the star.
!
! The simulation estimates the fraction of comets that:
!  (1) are ejected from the system
!  (2) collide with a target planet
!  (3) collide with a different planet
! Also estimates the mean cometary lifetime.
!
    module kinds
      integer, parameter::I4 = selected_int_kind(9)
      integer, parameter::R4 = kind(1.0)
      integer, parameter::R8 = kind(1.d0)
    end module kinds
!=========================================================
    module constants
      use kinds
      real(R8), parameter::ZERO  =  0.0_R8
      real(R8), parameter::ONE   =  1.0_R8
      real(R8), parameter::TWO   =  2.0_R8
      real(R8), parameter::THREE =  3.0_R8
      real(R8), parameter::FOUR  =  4.0_R8
      real(R8), parameter::FIVE  =  5.0_R8
      real(R8), parameter::SIX   =  6.0_R8
      real(R8), parameter::TEN   = 10.0_R8
!
      real(R8), parameter::HALF     = 0.5_R8
      real(R8), parameter::THIRD    = 1.0_R8 / 3.0_R8
      real(R8), parameter::TWOTHIRD = 2.0_R8 / 3.0_R8
!
      real(R8), parameter::PI     = 3.141592653589793_R8
      real(R8), parameter::TWOPI  = PI * TWO
      real(R8), parameter::PIBY2  = PI * HALF
      real(R8), parameter::ROOT2  = 1.41421356237309_R8
!
! Universal constants
      real(R8), parameter::GRAVCONST = 6.67428e-8_R8
!
! Solar System related constants
      real(R8), parameter::AU    = 1.49597870700e13_R8
      real(R8), parameter::GMSUN = 1.32712440041e26_R8
      real(R8), parameter::MSUN  = GMSUN / GRAVCONST
      real(R8), parameter::MEARTH   = MSUN / 332946.0487_R8
      real(R8), parameter::MJUPITER = MSUN / 1047.348644_R8
      real(R8), parameter::DAY  = 86400.0_R8
      real(R8), parameter::YEAR = DAY * 365.25_R8
!
! Maximum number of radial zones
      integer(I4), parameter::NMAX = 300000
!
! Typical comet eccentricity when apoapse is within
! the planetary system
      real(R8), parameter::ECC = 0.5_R8
!
! Distance of outermost bound orbit
      real(R8), parameter::ABOUND = 5e4_R8 * AU
!
! Target planet controls comet dynamics within NRH
! Hill radii of target planet
      real(R8), parameter::NRH = 5.0_R8
    end module constants
!=========================================================
    program markov
    use constants
    implicit none
    integer(I4)::i, j, k, n, itc_lo, itc_hi, itran
    integer(I4)::opt_target, opt_system
    real(R8)::a(NMAX), y(NMAX), torb(NMAX), mplan(NMAX)
    real(R8)::pcol(NMAX), pejc(NMAX), prip(NMAX)
    real(R8)::mat1(NMAX), mat2(NMAX), mat3(NMAX)
    real(R8)::rvec(NMAX), uvec(NMAX)
    real(R8)::mstar, ain, atarg, atc_lo, atc_hi, atran, aout
    real(R8)::mtarg, msmall, mtran, mout, density
    real(R8)::ftarg, fout, fcol1, fcol2, fejc1, fejc2, fejc3
    real(R8)::tbar, tmp
!---------------------------------------------------------
    call setup (n, k, a, y, torb, mplan, mstar, &
         ain, atarg, atc_lo, atc_hi, atran, aout, mtarg, &
         msmall, mtran, mout, density, opt_target, &
         opt_system, itc_lo, itc_hi, &
         itran, pcol, pejc, prip)
!
! Set up elements of a tridiagonal matrix
    mat1(1:n) = -HALF * (ONE - prip(1:n))
    mat2(1:n) = ONE
    mat3(1:n) = -HALF * (ONE - prip(1:n))
    mat2(1) = ONE - HALF * (ONE - prip(1))
    if (opt_target == 1) mat2(1) = ONE
!
! Fraction of comets that reach infinite distance
    rvec = ZERO
    rvec(n) = HALF * (ONE - prip(n))
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fout = uvec(k)
!
! Fraction of comets that hit target
    rvec = ZERO
    where (a < atc_hi.and.a > atc_lo) rvec = pcol
    if (opt_target == 1) rvec(1) = HALF * (ONE - prip(1))
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    ftarg = uvec(k)
!
! Fraction of comets that hit a small intermediate planet
    rvec = ZERO
    if (opt_system == 1) then
      where (a >= atc_hi.or.a < atc_lo) rvec = pcol
    else
      where (a < atran.and.(a >= atc_hi.or.a < atc_lo)) rvec = pcol
    end if
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fcol1 = uvec(k)
!
! Fraction of comets that hit a giant intermediate planet
    rvec = ZERO
    if (opt_system >= 2) then
      where (a >= atran) rvec = pcol
    end if
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fcol2 = uvec(k)
!
! Fraction of comets ejected by target
    rvec = ZERO
    where (a < atc_hi.and.a > atc_lo) rvec = pejc
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fejc1 = uvec(k)
!
! Comets with a < a_out ejected by an intermediate planet
    rvec = ZERO
    where (a < aout.and.(a >= atc_hi.or.a < atc_lo)) rvec = pejc
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fejc2 = uvec(k)
!
! Comets with a > a_out ejected by an intermediate planet
    rvec = ZERO
    where (a >= aout) rvec = pejc
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    fejc3 = uvec(k)
!
! Tranpose the tridiagonal matrix
    mat1(2:n)   = -HALF * (ONE - prip(1:n-1))
    mat2(1:n)   = ONE
    mat3(1:n-1) = -HALF * (ONE - prip(2:n))
    mat2(1) = ONE - HALF * (ONE - prip(1))
    if (opt_target == 1) mat2(1) = ONE
!
! Calculate the mean cometary lifetime
! (sum of kth row of inverted matrix weighted by TORB)
    rvec = ZERO
    rvec(k) = ONE
    call tridag (mat1, mat2, mat3, rvec, uvec, n)
    tbar = sum (uvec(1:n) * torb(1:n))
!
! Output the results
    call output (mstar, ain, atarg, atc_lo, atc_hi, &
         atran, aout, mtarg, msmall, mtran, mout, density, &
         opt_target, opt_system, n, k, itc_lo, itc_hi, &
         itran, fejc1, fejc2, &
         fejc3, fout, ftarg, fcol1, fcol2, tbar)
!
    end program markov
!=========================================================
    subroutine calc_pcol_pejc (n, a, y, mplan, mstar, &
         aout, density, pcol, pejc, prip)
    use constants
    implicit none
    integer(I4), intent(in)::n
    real(R8), intent(in)::y(n), a(n), mplan(n)
    real(R8), intent(in)::mstar, aout, density
    real(R8), intent(out)::pcol(n), pejc(n), prip(n)
!
    integer(I4)::i
    real(R8)::ap(n), rp(n), mscl(n), e2(n), e3(n)
!---------------------------------------------------------
    ap = min (a, aout)
    rp = (THREE * mplan / (FOUR * PI * density))**THIRD
    mscl = mplan / mstar
    e2 = ECC * ECC
    e3 = ECC * e2
!
    pcol = rp * rp / (ap * ap * ECC) &
          * (ONE  +  TWO * mscl * ap / (rp * e2))
    pejc = FOUR * mscl * mscl * a * a / (ap * ap * e3)
!
! Make sure probabilities don't exceed unity
    pcol = min (pcol, ONE)
    pejc = min (pejc, ONE)
!
! Note: collision and ejection cross sections overlap
    prip = max (pcol, pejc)
    pejc = prip - pcol
!
    end subroutine calc_pcol_pejc
!=========================================================
    subroutine output (mstar, ain, atarg, atc_lo, atc_hi, &
         atran, aout, mtarg, msmall, mtran, mout, density, &
         opt_target, opt_system, n, k, itc_lo, itc_hi, &
         itran, fejc1, fejc2, &
         fejc3, fout, ftarg, fcol1, fcol2, tbar)
    use constants
    implicit none
    integer(I4), intent(in)::n, k, itc_lo, itc_hi, itran
    integer(I4), intent(in)::opt_target, opt_system
    real(R8), intent(in)::mstar, ain, atarg, atc_lo, atc_hi
    real(R8), intent(in)::atran, aout, mtarg
    real(R8), intent(in)::msmall, mtran, mout
    real(R8), intent(in)::density, tbar
    real(R8), intent(in)::fejc1, fejc2, fejc3, fout
    real(R8), intent(in)::ftarg, fcol1, fcol2
!
    integer(I4)::i, j
    character(26)::c
!---------------------------------------------------------
    open (21, file='log.out', status='replace')
    do i = 1, 2
      j = 6
      if (i == 2) j = 21
!
      write (j,*)
      write (j,203) 'Stellar mass (solar):             ', &
           mstar / MSUN
      write (j,201) 'Innermost planet distance (AU):   ', &
           ain / AU
      write (j,201) 'Distance of target (AU):          ', &
           atarg / AU
      write (j,201) 'Inner edge of target control (AU):', &
           atc_lo / AU
      write (j,201) 'Outer edge of target control (AU):', &
           atc_hi / AU
      write (j,201) 'Small/giant planet transition(AU):', &
           atran / AU
      write (j,201) 'Outermost giant-plan dist (AU):   ', &
           aout / AU
      write (j,*)
      write (j,201) 'Target planet mass (Earth):       ', &
           mtarg / MEARTH
      write (j,201) 'Small-planet mass (Earth):        ', &
           msmall / MEARTH
      write (j,201) 'Innermost giant-plan mass (Earth):', &
           mtran / MEARTH
      write (j,201) 'Outermost giant-plan mass (Earth):', &
           mout / MEARTH
      write (j,201) 'Planetary density (g/cm^3):       ', &
           density
!
      c = 'Target: ordinary planet   '
      if (opt_target == 1) c = 'Target: absorbing boundary'
      write (j,'(1x,a)') c
!
      c = 'giant and small planets   '
      if (opt_system == 1) c = 'small planets only        '
      if (opt_system == 2) c = 'giant planets only        '
      write (j,'(1x,2a)') 'System: ', c
!
201   format (1x,a,1x,es10.3)
203   format (1x,a,f7.3)
!
      write (j,*)
      write (j,204) 'Number of radial zones:          ', n
      write (j,204) 'Starting zone for comets:        ', k
      write (j,204) 'Small/giant plan transition zone:', itran
      write (j,204) 'Inner edge of target control:    ', itc_lo
      write (j,204) 'Outer edge of target control:    ', itc_hi
204   format (1x,a,1x,i6)
!
      write (j,*)
      write (j,205) 'Comets that hit the target:       ', &
           ftarg
      write (j,205) 'Hit another small planet:         ', &
           fcol1
      write (j,205) 'Hit another giant plannet:        ', &
           fcol2
      write (j,205) 'Ejected by the target:            ', &
           fejc1
      write (j,205) 'Otherwise ejected when a < a_out: ', &
           fejc2
      write (j,205) 'Otherwise ejected when a > a_out: ', &
           fejc3
      write (j,205) 'Reached outer edge of the grid:   ', &
           fout
      write (j,205) 'Total probability:                ', &
           ftarg + fcol1 + fcol2 + fejc1 + fejc2 + fejc3 + fout
205   format(1x,a,1x,es10.3)
!
      write (j,*)
      write (j,206) 'Mean comet lifetime (year):', tbar / YEAR
206   format(1x,a,1x,es10.3)
      write (j,*)
    end do
!
    close (21)
!
    end subroutine output
!=========================================================
    subroutine setup (n, k, a, y, torb, mplan, mstar, &
         ain, atarg, atc_lo, atc_hi, atran, aout, mtarg, &
         msmall, mtran, mout, density, opt_target, &
         opt_system, itc_lo, itc_hi, &
         itran, pcol, pejc, prip)
    use constants
    implicit none
    integer(I4), intent(out)::n, k, itc_lo, itc_hi, itran
    integer(I4)::opt_target, opt_system
    real(R8), intent(out)::a(NMAX), y(NMAX), torb(NMAX)
    real(R8), intent(out)::mplan(NMAX)
    real(R8), intent(out)::ain, atarg, atran, aout
    real(R8), intent(out)::msmall, mtran, mout, atc_lo
    real(R8), intent(out)::mstar, density, mtarg, atc_hi
    real(R8), intent(out)::pcol(NMAX),pejc(NMAX), prip(NMAX)
!
    integer(I4)::i
!---------------------------------------------------------
    open (20, file='markov.in', status='old')
    read (20,*) mstar
    read (20,*) ain
    read (20,*) atarg
    read (20,*) atran
    read (20,*) aout
!
    read (20,*) mtarg
    read (20,*) msmall
    read (20,*) mtran
    read (20,*) mout
!
    read (20,*) density
    read (20,*) opt_target
    read (20,*) opt_system
    close (20)
!
! Convert to cgs units
    mstar  = mstar * MSUN
    mtarg  = mtarg * MEARTH
    msmall = msmall * MEARTH
    mtran  = mtran * MEARTH
    mout   = mout  * MEARTH
    ain   = ain   * AU
    atarg = atarg * AU
    atran = atran * AU
    aout  = aout  * AU
!
    call setup_radial_grid (mstar, ain, atarg, atc_lo, &
         atc_hi, atran, aout, mtarg, msmall, mtran, mout, &
         itc_lo, itc_hi, itran, n, k, a, y, torb, mplan, &
         opt_target, opt_system)
!
! Get probability of collision and ejection by planets
    call calc_pcol_pejc (n, a, y, mplan, mstar, &
         aout, density, pcol, pejc, prip)
!
    open (20, file='zones.out', status='replace')
    write (20,*) '    i        y          a (AU)  ', &
         '   P (y)    Mplan (M_E)     Pcol        ', &
         'Pejc   '
    write (20,*) '--------------------------------', &
         '----------------------------------------', &
         '-------'
    do i = 1, n
      write(20,201) i, y(i), a(i) / AU, &
           torb(i) / YEAR, mplan(i) / MEARTH, &
           pcol(i), pejc(i)
    end do
    close (20)
201 format (1x,i6,1x,f10.5,1x,f13.5,4(1x,es11.4))
!
    end subroutine setup
!=========================================================
    subroutine setup_radial_grid (mstar, ain, atarg, atc_lo, &
         atc_hi, atran, aout, mtarg, msmall, mtran, mout, &
         itc_lo, itc_hi, itran, n, k, a, y, torb, mplan, &
         opt_target, opt_system)
    use constants
    implicit none
    integer(I4), intent(in)::opt_target, opt_system
    real(R8), intent(in)::mstar, ain, atarg, aout
    real(R8), intent(in)::mtarg, msmall, mtran
    real(R8), intent(inout)::atran, mout
    integer(I4), intent(out)::n, k, itc_lo, itc_hi, itran
    real(R8), intent(out)::a(NMAX), y(NMAX), torb(NMAX)
    real(R8), intent(out)::mplan(NMAX), atc_lo, atc_hi
!
    integer(I4)::i
    real(R8)::h, expo, yin, ytran, ytc_lo, ytc_hi
    real(R8)::ytarg, ybound, dy, tn
!---------------------------------------------------------
! Target planet controls region within NRH Hill radii
    h = (mtarg * THIRD / mstar)**THIRD
    atc_lo = atarg * (ONE  -  NRH * h)
    atc_hi = atarg * (ONE  +  NRH * h)
    atc_lo = max (atc_lo, ain)
    atran = max (atran, atc_hi)
!
    if (opt_target == 1) atc_hi = atarg
!
! If there are no giant planets
    if (opt_system == 1) then
      atran = aout
      expo = ZERO
      mout = msmall
    else
      expo = log(mtran / mout) / log(atran / aout)
    end if
!
! If there are no small planets
    if (opt_system == 2) atran = atc_hi
!
! The Y values of some important states
    if (expo == ZERO) then
      ytran = log(atran / aout)
    else
      ytran = (ONE - (aout / atran)**expo) / expo
    end if
!
    ytc_hi = ytran  +  mout / msmall * log(atc_hi / atran)
    ytarg  = ytc_hi +  mout / mtarg  * log(atarg / atc_hi)
    ytc_lo = ytarg  +  mout / mtarg  * log(atc_lo / atarg)
    yin    = ytc_lo +  mout / msmall * log(ain / atc_lo)
    ybound = ONE  -  aout / ABOUND
!
! Width of each non-abosrbing state in terms of Y
    dy = TEN * mout / mstar
!
! Number of non-absorbing states
    n = (ybound - yin) / dy
!
! Starting state for comets
    k = -yin / dy
!
! Y values and distances of each non-absorbing state
    do i = 1, n
      y(i) = yin  +  dy * i
!
      if (y(i) >= ZERO) then
        a(i) = aout / (ONE - y(i))
      else if (y(i) < ytc_lo) then
        a(i) = atc_lo * exp(msmall * (y(i) - ytc_lo) / mout)
      else if (y(i) < ytc_hi) then
        a(i) = atc_hi * exp(mtarg  * (y(i) - ytc_hi) / mout)
      else if (y(i) < ytran) then
        a(i) = atran * exp(msmall * (y(i) - ytran) / mout)
      else
        if (expo == ZERO) then
          a(i) = aout * exp(y(i))
        else
          a(i) = aout / (ONE  -  expo * y(i))**(ONE/expo)
        end if
      end if
    end do
!
! Indices of some important states
    do i = n, 1, -1
      if (atran < a(i)) itran = i
      if (atc_hi < a(i)) itc_hi = i
      if (atc_lo < a(i)) itc_lo = i
    end do

    write (*,*) 'i, ytran: ', itran, ytran
    write (*,*) 'i, ytc_hi ', itc_hi, ytc_hi
    write (*,*) 'i, ytc_lo ', itc_lo, ytc_lo
    write (*,*) 'i, yin    ', 1, yin
!
! Orbital period of outermost planet
    tn = TWO * PI * sqrt(aout**3 / GRAVCONST / mstar)
!
! Orbital period for each state
    torb = tn * (a / aout)**1.5d0
!
! The mass of the controlling planet for each state
    where (a > aout)
      mplan = mout
    else where (atc_lo < a.and. a < atc_hi)
      mplan = mtarg
    else where (a < atran)
      mplan = msmall
    else where
      mplan = mout * (a / aout)**expo
    end where
!
    end subroutine setup_radial_grid
!=========================================================
! Solves for a vector U where R = T U, where R is a
! vector and T is a tridiagonal matrix with elements
! given in A, B and C.
! Based on a routine from Numerical Recipes, Press et al.
    subroutine tridag (a, b, c, r, u, n)
    use constants
    implicit none
    integer(I4), intent(in)::n
    real(R8), intent(in)::a(n), b(n), c(n), r(n)
    real(R8), intent(out)::u(n)
!
    integer(I4)::j
    real(R8)::bet, gam(NMAX)
!---------------------------------------------------------
    bet = b(1)
    u(1) = r(1) / bet
!
    do j = 2, n
      gam(j) = c(j-1) / bet
      bet = b(j) - a(j) * gam(j)
      u(j) = (r(j) - a(j) * u(j-1)) / bet
    end do
!
    do j = n - 1, 1, -1
      u(j) = u(j) - gam(j+1) * u(j+1)
    end do
!
    end subroutine tridag



 
    
