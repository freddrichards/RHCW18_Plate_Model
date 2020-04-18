program RHCW18

! plate cooling model with pressure- and temperature-dependent thermal
! parameters and an insulating crustal layer from the paper:
!
! Fred Richards, Mark Hoggard, Laurence Cowton & Nicky White (2018), 
! Reassessing the thermal structure of oceanic lithosphere with revised 
! global inventories of basement depths and heat flow measurements,
! Journal of Geophysical Research: Solid Earth, Volume 123, 9136-9161.
!
! v1.1 - August 2019
  
  use global
  implicit none
  
  ! ------------------------------------------------------------------------
  ! Set input parameters
  ! ------------------------------------------------------------------------
  
  ! model setup
  integer,  parameter :: itmax = 1e8          ! maximum number of iterations
  real(wp), parameter :: secMa = 3.1556736e13 ! seconds per Ma
  real(wp), parameter :: agemax  = 300 		  ! model run time in Ma
  real(wp), parameter :: dzinit = 1000.       ! initial depth interval between nodes (i and i+1 positions)
  integer,  parameter :: zc     = 7           ! integer - node corresponding to crustal thickness (must be nearest integer number of depth node spacings)
  real(wp), parameter :: T0=273.              ! reference temperature
  real(wp), parameter :: roundtol=1.          ! rounding tolerance (K)
  real(wp), parameter :: dageplot = 0.1       ! time between outputs (Ma)



  ! reference densities
  real(wp), parameter :: p0=3330. ! mantle
  real(wp), parameter :: p0c=2950. ! crust
  real(wp), parameter :: pw0=1030. ! water
  real(wp), parameter :: g=9.81 ! gravitational acceleration
  
  ! pressure-dependent density (Grose & Afonso, 2013)
  real(wp), parameter ::  K0=130.e9 ! bulk modulus at reference temperature and pressure
  real(wp), parameter ::  KT=4.8 ! temperature derivative of bulk modulus
  real(wp), parameter ::  grun=6. ! Gruneisen parameter
  
  
  
  ! olivine thermal expansivity (Bouhifd et al. 1996)
  real(wp), parameter :: a0=2.832e-5
  real(wp), parameter :: a1=0.758e-8
  
  ! crustal thermal expansivity Grose & Afonso, 2013)
  ! diopside
  real(wp), parameter :: a0di=2.1e-5 
  real(wp), parameter :: a1di=1.75e-8
  ! anorthite
  real(wp), parameter :: a0an=0.87e-5
  real(wp), parameter :: a1an=0.9e-8
  ! albite
  real(wp), parameter :: a0ab=1.75e-5
  real(wp), parameter :: a1ab=1.95e-8
  ! crust
  real(wp), parameter :: a0c=1.639e-5
  real(wp), parameter :: a1c=1.322e-8



  ! mantle heat capacity (Korenaga & Korenaga 2016): 
  real(wp), parameter :: c0_mantle=1580.
  real(wp), parameter :: c1_mantle =12230.
  real(wp), parameter :: c2_mantle =1694.e6
  
  ! crustal heat capacity (Grose and Afonso, 2013)
  ! olivine
  real(wp), parameter :: k0fofa=1.6108e3
  real(wp), parameter :: k1fofa=-1.24788e4
  real(wp), parameter :: k3fofa=-1.728477e9
  ! clinopyroxene
  real(wp), parameter :: c0_cpx=2.1715e3
  real(wp), parameter :: c1_cpx=-4.555e-1
  real(wp), parameter :: c2_cpx=1.1332e6
  real(wp), parameter :: c3_cpx=-2.22716e4
  real(wp), parameter :: c4_cpx=1.299e-4
  ! plagioclase
  real(wp), parameter :: c0_plag =1.85757e3
  real(wp), parameter :: c1_plag =-3.324e-1
  real(wp), parameter :: c2_plag =-5.061e6
  real(wp), parameter :: c3_plag =-1.64946e4
  real(wp), parameter :: c4_plag =1.505e-4




  ! mantle diffusivity (Grose & Afonso, 2013)
  real(wp), parameter ::  a_ol =0.565
  real(wp), parameter ::  b_ol =0.67
  real(wp), parameter ::  c_ol =590.
  real(wp), parameter ::  d_ol =1.4
  real(wp), parameter ::  e_ol =135
 
  ! crustal diffusivity (Grose & Afonso, 2013)
  ! clinopyroxene
  real(wp), parameter ::  a_cpx =0.59
  real(wp), parameter ::  b_cpx =1.03 
  real(wp), parameter ::  c_cpx =386. 
  real(wp), parameter ::  d_cpx =0.928
  real(wp), parameter ::  e_cpx =125.
  ! plagioclase
  real(wp), parameter ::  a_plag =0.36
  real(wp), parameter ::  b_plag =0.4
  real(wp), parameter ::  c_plag =300.
  !    crust
  real(wp), parameter ::  a_crust =0.414
  real(wp), parameter ::  b_crust =0.22
  real(wp), parameter ::  c_crust =455.
  real(wp), parameter ::  d_crust =0.371
  real(wp), parameter ::  e_crust =226.
  
  ! Pressure dependence of lattice conductiviity
  real(wp), parameter ::  dkP =0.05 ! in GPa^-1
  
  ! Radiative conductivity (Grose & Afonso, 2013)
  real(wp), parameter ::  d = 0.5
  real(wp), parameter ::  Ar = (1.8*(1.-exp(-((d**(1.3))/0.15))))-(1.-exp(-((d**(0.5))/5.)))
  real(wp), parameter ::  Br = (11.7*exp(-d/0.159))+(6.*exp(-((d**(3.))/10.)))
  real(wp), parameter ::  Tar = 490.+(1850.*exp(-((d**0.315)/0.825)))+(875.*exp(-d/0.18))
  real(wp), parameter ::  Tbr = 2700.+(9000.*exp(-((d**0.5)/0.205)))
  real(wp), parameter ::  xa = 167.5+(505.*exp(-((d**0.5)/0.85)))
  real(wp), parameter ::  xb = 465.+(1700.*exp(-((d**0.94)/0.175))) 
  
  ! ------------------------------------------------------------------------
  ! Read input arguments
  ! ------------------------------------------------------------------------
  integer :: narg, cptArg, zp, rd, nz, Tpot ! no. of arguments, argument counter and input arguments
  character(30)  ::  name
  
  narg=command_argument_count() ! count number of input arguments
  !Loop over the arguments
  if(narg==3)then
    do cptArg=1,narg
      call get_command_argument(cptArg,name)
      if(cptArg==1)then
        read(name,*)Tpot ! read potential temperature in deg C
      else if(cptArg==2)then
      read(name,*)zp ! read plate thickness in km
      else if(cptArg==3)then
      read(name,*)rd ! read ridge depth in m
      end if
    end do
  else
    ! Throw error if wrong number of input arguments.
    write(*,*)"Incorrect number of input parameters, format is ./cool [Tp] [zp] [zr]"
    stop
  end if 
  
  nz=zp ! set plate thickness to number of data points
  
  call plate_cooling(Tpot,nz) ! call main program
  
  
! ------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------  
  
  subroutine plate_cooling(Tpot,nz)
    ! ------------------------------------------------------------------------
    ! Master subroutine
    ! ------------------------------------------------------------------------
    use global
    implicit none

    ! set variables
    real(wp) :: t, dt, Tbase, dtinit, tplot
    real(wp) :: sump0, ref_rd, ref_sub, xmin
    integer :: it, i, nz, Tpot
    real(wp), dimension(0:nz)  :: Tm, Told, pdat0, prdat0, dz0
    real(wp), dimension(0:nz)  :: dV0, alphaP0, rhoP0, intalphaT0, z0, misfit
    
    ! initialise filenames
    character(30)  ::  hffile
    character(30)  ::  Tinfile
    character(30)  ::  Tgridfile
    character(30)  ::  depthfile
    
    ! set file names
    if (nz<100) then
      write (hffile, "('hf-',I4,'-',I2,'.dat')")Tpot,zp ! output file for heat flow
      write (depthfile, "('depth-',I4,'-',I2,'-',I4,'.dat')")Tpot,zp,rd ! output file for subsidence
      write (Tinfile, "('mac',I4,'.dat')")Tpot ! initial temperature profile
      write (Tgridfile, "('Tgrid-',I4,'-',I2,'.dat')")Tpot,zp ! output file for thermal structure
    else 
      write (hffile, "('hf-',I4,'-',I3,'.dat')")Tpot,zp ! output file for heat flow
      write (depthfile, "('depth-',I4,'-',I3,'-',I4,'.dat')")Tpot,zp,rd ! output file for subsidence
      write (Tinfile, "('mac',I4,'.dat')")Tpot ! initial temperature profile
      write (Tgridfile, "('Tgrid-',I4,'-',I3,'.dat')")Tpot,zp ! output file for thermal structure
    end if
    
    ! open initial temperature profile
    open(4, file = Tinfile)

    ! read in intital temperature profile
    do i=0,nz
      read(4,*) Tm(i)            ! melting parameterisation-derived temperature in deg C
    end do
    Tm=Tm+273.                   ! convert to K
    close(4)
    
    ! set initial values
    Tbase=Tm(nz)                ! set basal temperature to adiabatic value at zp
    t = 0._wp                   ! set time = 0
	tplot = dageplot*secMa	    ! set first datapoint in time
    dtinit=(dzinit**2.)/(2.2*(6./((3100.*100.)/0.14069)))  ! set timestep (see McKenzie et al., 2005)
    dt = dtinit ! set timestep to initial value
    Told=Tm ! set initial temperature profile
    ref_rd=rd/1000. ! reference ridge depth used in initial pressure profile calculation
    ref_sub=0. ! reference subsidence at t=0
    
    ! calculate initial pressure profile
    ! crust
    do i=0,zc
      dz0(i)=dzinit ! initial depth-spacing of nodes
      z0(i)=i*dz0(i) ! initial. depth within plate
      if (i.eq.0) then
          prdat0(i)=(1028.+(2.4*(ref_sub+ref_rd)))*g*ref_rd*dz0(i) ! pressure at top of plate from water column
      else
          prdat0(i)=prdat0(i-1)+((pdat0(i-1)*g)*dz0(i)) ! pressure at successive nodes 
      endif    
      misfit(i)=brent(xmin,prdat0(i)) ! determine (V0/V)T
      dV0(i)=xmin ! set (V0/V)T
      alphaP0(i)=dV0(i)*exp((grun+1)*((dV0(i)**(-1.))-1.)) ! calculate pressure-dependent expansivity
      rhoP0(i)=p0c*dV0(i) ! calculate pressure dependent density
      intalphaT0(i)=(a0c*(Tm(i)-T0))+((a1c/2.)*((Tm(i)**2.)-(T0**2.))) ! calculate thermal expansivity
      pdat0(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT0(i))) ! calculate pressure and temperature-dependent density
    end do
    ! mantle
    do i=zc+1,nz
      dz0(i)=dzinit ! initial depth-spacing of nodes
      z0(i)=i*dz0(i) ! initial. depth within plate
      prdat0(i)=prdat0(i-1)+((pdat0(i-1)*g)*dz0(i)) ! pressure at successive nodes 
      misfit(i)=brent(xmin,prdat0(i)) ! determine (V0/V)T
      dV0(i)=xmin ! set (V0/V)T
      alphaP0(i)=dV0(i)*exp((grun+1)*((dV0(i)**(-1.))-1.)) ! calculate pressure-dependent expansivity
      rhoP0(i)=p0*dV0(i) ! calculate pressure dependent density
      intalphaT0(i)=(a0*(Tm(i)-T0))+((a1/2.)*((Tm(i)**2.)-(T0**2.))) ! calculate thermal expansivity
      pdat0(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT0(i))) ! calculate pressure and temperature-dependent density
    end do
    sump0=sum(pdat0) ! total density
    
    ! open files for writing output
    open(2, file = hffile, position = 'append')  
    open(5, file = Tgridfile, position = 'append')
    open(7, file = depthfile, position = 'append')

    ! output data for t=0
    call output(Tm,t,sump0,alphaP0,rhoP0,pdat0,prdat0)
    
    ! set iteration counter
    it = 1

    ! enter time loop
    do while ((t/secMa < agemax) .and. (it < itmax)) ! while within timeframe and max iterations
  	  call timestep(Tm,Told,t,dt,Tbase,alphaP0,rhoP0,prdat0,pdat0,it) ! finite difference calculation step
      t = t+dt ! reset time
      if (t > tplot) then		! if output required
 	     call output(Tm,t,sump0,alphaP0,rhoP0,pdat0,prdat0) ! write output to file
  	     tplot = tplot + (dageplot*secMa)			        ! increase tplot by dtplot
      end if
      it = it + 1                            ! add 1 to iteration count
    end do
    
    ! close output files after writing
    close(2)
    close(5)
    close(7)

  end subroutine plate_cooling
  
! ------------------------------------------------------------------------
    
  subroutine  output(Ta,t,sump0,alphaP0,rhoP0,pdat0,prdat0)
    ! ------------------------------------------------------------------------
    ! Model output subroutine
    ! ------------------------------------------------------------------------
    use global
    implicit none

    real(wp), dimension(0:nz), intent(in) :: Ta,pdat0,alphaP0,rhoP0,prdat0
    real(wp), intent(in) :: t,sump0
    real(wp), dimension(0:nz) :: intalphaT,pdat,cpdat,Ddat,krad,cum_contract,mean_contract,dz,cumz
    real(wp) :: age,ktot,H,sump,depth,xmin,contraction,pb,comp,j,pw,rdk
    integer  :: i,minflag

    ! convert time to seconds
    age=t/secMa

    minflag=0 ! set flag for working out minimum depth at which densities match between timesteps
   
    ! calculate isostatic depth
    comp=0. ! set initial compensation depth to zero
    pb=0. ! set initial density at compensation depth to zero
    ! crust
    do i=0,zc
      intalphaT(i)=(a0c*(Ta(i)-T0))+((a1c/2.)*((Ta(i)**2.)-(T0**2.))) ! calculate temperature dependence of expansivity
      pdat(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT(i))) ! calculate pressure and temperature-dependent density
      if (i.eq.0) then
        cumz(i)=0. ! set cumulative depth to zero
        dz(i)=0. ! set initial depth-spacing of nodes to zero
      else
        dz(i)=(((pdat0(i-1)/pdat(i-1))+(pdat0(i)/pdat(i)))/2.)*dzinit ! update depth spacing to account for vertical contraction (assumes vertical contraction so pressure at base of each cell is conserved)
          cumz(i)=cumz(i-1)+dz(i) ! update cumulative depth
      endif
      ! if at minimum depth where densities from successive timesteps match to nearest integer, define compensation depth and corresponding density
      if (i.gt.0.and.(int(pdat0(i)*roundtol)/roundtol).eq.(int(pdat(i)*roundtol)/roundtol).and.t.gt.0.and.minflag.eq.0) then
          comp=cumz(i) ! set compensation depth
          pb=pdat(i) ! set density at compensation depth
          minflag=minflag+1 ! reset minimum depth flag
      endif
    end do
    ! mantle
    do i=zc+1,nz
      intalphaT(i)=(a0*(Ta(i)-T0))+((a1/2.)*((Ta(i)**2.)-(T0**2.))) ! calculate temperature dependence of expansivity
      pdat(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT(i))) ! calculate pressure and temperature-dependent density
      dz(i)=(((pdat0(i-1)/pdat(i-1))+(pdat0(i)/pdat(i)))/2.)*dzinit ! update depth spacing to account for vertical contraction (assumes vertical contraction so pressure at base of each cell is conserved)
      cumz(i)=cumz(i-1)+dz(i) ! update cumulative depth
      ! if at minimum depth where densities from successive timesteps match to nearest integer, define compensation depth and corresponding density
      if ((int(pdat0(i)*roundtol)/roundtol).eq.(int(pdat(i)*roundtol)/roundtol).and.t.gt.0.and.minflag.eq.0) then
          comp=cumz(i) ! set compensation depth
          pb=pdat(i) ! set density at compensation depth
          minflag=minflag+1 ! reset minimum depth flag
      endif
    end do

    ! determine total contraction of column
    sump=sum(pdat) ! cumulative density
    contraction=(1.-(sump0/sump))*zp*dzinit ! total contraction (assumes vertical contraction so pressure at base of each cell is conserved)
 
    ! determine subsidence as a function of age using density at compensation depth and pressure-dependence of water density
    rdk=rd/1000. ! ridge depth in km
    depth=(rtbis(pb,(contraction/dzinit),pw,rdk)+rdk)*dzinit !iterates to find subsidence as water density is a function of depth 
    write(7,'(2e25.16e3)') age, depth ! write out age (Ma) and depth (m)
    
    ! determine temperature evolution as a function of age and depth within the plate
    cum_contract(0)=1.
    j=1.
    do i=0,nz
      cum_contract(i)=((pdat0(i)/pdat(i))+cum_contract(i-1)) ! cumulative contraction (assumes vertical contraction so pressure at base of each cell is conserved)
      mean_contract(i)=cum_contract(i)/j ! mean contraction
      write(5,'(4e25.16e3)') age, i*mean_contract(i), (i*mean_contract(i))+(depth/dzinit), Ta(i)-T0 ! write out age (Ma), depth below top of plate (km), depth below sea-level (km) and temperature (°C).
      j=j+1.
    end do

    ! determine heat flow as a function of age
    cpdat(0)=((0.15*((k0fofa)+(k1fofa*(Ta(0)**(-1./2.)))+(k3fofa*(Ta(0)**(-3.)))))+&
		(0.2*(c0_cpx+(c1_cpx*Ta(0))+(c2_cpx*(Ta(0)**(-2)))+(c3_cpx*(Ta(0)**(-1./2.)))+&
			(c4_cpx*(Ta(0)**2.))))+(0.65*(c0_plag+(c1_plag*Ta(0))+(c2_plag*(Ta(0)**(-2)))+&
				(c3_plag*(Ta(0)**(-1./2.)))+(c4_plag*(Ta(0)**2.))))) ! heat capacity 
    Ddat(0)=(a_crust+(b_crust*exp(-(Ta(0)-T0)/c_crust))+(d_crust*exp(-(Ta(0)-T0)/e_crust)))*1e-6 ! diffusivity
    krad(0)=Ar*exp(-(((Ta(0)-Tar)**(2.))/(2*(xa**(2.)))))+Br*exp(-(((Ta(0)-Tbr)**(2.))/(2*(xb**(2.))))) ! radiative conductivity
    ktot=(pdat(0)*Ddat(0)*cpdat(0)*exp(dkP*(prdat0(0)/1.e+9)))+krad(0) ! total conductivity
    H=(ktot*(Ta(1)-Ta(0)))/dz(1) ! heat flow
    write(2,'(2e25.16e3)') age, H

  end subroutine output
  
! ------------------------------------------------------------------------
  
  subroutine timestep(Tt1,Tt3,t,twodt,Tbase,alphaP0,rhoP0,prdat0,pdat0,it)

    ! predictor-corrector routine (n is timestep, i is position)
    
    use global
    implicit none
    real(wp), dimension(0:nz), intent(inout) :: Tt1, Tt3, alphaP0,rhoP0,prdat0,pdat0
    real(wp), intent(in) :: t, twodt, Tbase
    integer, intent(in) :: it
    real(wp), dimension(0:nz) :: Tt2, cz, czm, cora, corb, p, cp, dz

    ! predicts the T at n+1 using knowledge of T at n
    ! then corrects estimate of k(n+1/2) using average of T(n+1) and T(n) and recalculates T at n+1
    ! [A][x] = [b], where matrix A relates to prefactors of T(n+1,:), and b to T(n,:)
    ! This routine is not worth doing more than once per timestep

    ! predictor step    
    ! Tt3 = n-1 timestep, Tt1 = n timestep (present), Tt2 = n+1 timestep (initial guess future)    
    call correction(Tt3,Tt1,Tt1,cora,alphaP0,rhoP0) ! optional dpCp/dt correction term 
    call set_var(Tt1,cz,czm,p,cp,alphaP0,rhoP0,prdat0,pdat0,dz) ! set variables
    call set_eqn(Tt1,Tt2,cz,czm,p,cp,cora,twodt,Tbase,dz) ! solve [A][x] = [b]

    ! set Tt3 to Tt1 so it n timestep becomes the n-1 timestep next time the routine is called
    Tt3=Tt1

    ! corrector step
    ! at the end of this step The n+1 timestep result becomes Tt1 (second variable of the set_eqn call is the output) 
    call correction(Tt1,Tt2,Tt1,corb,alphaP0,rhoP0) ! optional dpCp/dt correction term 
    call set_var((Tt1+Tt2)/2.,cz,czm,p,cp,alphaP0,rhoP0,prdat0,pdat0,dz) ! set variables
    call set_eqn(Tt1,Tt1,cz,czm,p,cp,corb,twodt,Tbase,dz) ! solve [A][x] = [b]


  end subroutine timestep
  
! ------------------------------------------------------------------------
  
subroutine correction(T1,T2,T3,cora,alphaP0,rhoP0)
    
    ! finds  optional dpCp/dt correction term used required in McKenzie et al. 2005 eqn. 3

    use global
    implicit none
    real(wp), dimension(0:nz), intent(in)  ::  T1,T2,T3,alphaP0,rhoP0
    real(wp), dimension(0:nz), intent(out)  :: cora
    real(wp), dimension(0:nz) :: p1,cp1,p2,cp2,p3,cp3
    real(wp), dimension(0:nz) :: intalphaT1,intalphaT2,intalphaT3
    integer  ::  i
    
    do i=0,zc
      intalphaT1(i)=(a0c*(T1(i)-T0))+((a1c/2.)*((T1(i)**2.)-(T0**2.)))
      p1(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT1(i)))
      cp1(i)=((0.15*((k0fofa)+(k1fofa*(T1(i)**(-1./2.)))+(k3fofa*(T1(i)**(-3.)))))+&
	  	(0.2*(c0_cpx+(c1_cpx*T1(i))+(c2_cpx*(T1(i)**(-2.)))+(c3_cpx*(T1(i)**(-1./2.)))+&
			(c4_cpx*(T1(i)**2.))))+(0.65*(c0_plag+(c1_plag*T1(i))+(c2_plag*(T1(i)**(-2.)))+&
				(c3_plag*(T1(i)**(-1./2.)))+(c4_plag*(T1(i)**2.)))))
    end do
    do i=zc+1,nz
      intalphaT1(i)=(a0*(T1(i)-T0))+((a1/2.)*((T1(i)**2.)-(T0**2.)))
      p1(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT1(i)))
      cp1(i)=(c0_mantle-(c1_mantle*(T1(i)**(-1./2.)))-(c2_mantle*(T1(i)**(-3.))))
    end do

    do i=0,zc
      intalphaT2(i)=(a0c*(T2(i)-T0))+((a1c/2.)*((T2(i)**2.)-(T0**2.)))
      p2(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT2(i)))
      cp2(i)=((0.15*((k0fofa)+(k1fofa*(T2(i)**(-1./2.)))+(k3fofa*(T2(i)**(-3.)))))+&
      	(0.2*(c0_cpx+(c1_cpx*T2(i))+(c2_cpx*(T2(i)**(-2.)))+(c3_cpx*(T2(i)**(-1./2.)))+&
      		(c4_cpx*(T2(i)**2.))))+(0.65*(c0_plag+(c1_plag*T2(i))+(c2_plag*(T2(i)**(-2.)))+&
        		(c3_plag*(T2(i)**(-1./2.)))+(c4_plag*(T2(i)**2.)))))
    end do
    do i=zc+1,nz
      intalphaT2(i)=(a0*(T2(i)-T0))+((a1/2.)*((T2(i)**2.)-(T0**2.)))
      p2(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT2(i)))
      cp2(i)=(c0_mantle-(c1_mantle*(T2(i)**(-1./2.)))-(c2_mantle*(T2(i)**(-3.))))
    end do
    
    do i=0,zc
      intalphaT3(i)=(a0c*(T3(i)-T0))+((a1c/2.)*((T3(i)**2.)-(T0**2.)))
      p3(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT3(i)))
      cp3(i)=((0.15*((k0fofa)+(k1fofa*(T3(i)**(-1./2.)))+(k3fofa*(T3(i)**(-3.)))))+&
	  	(0.2*(c0_cpx+(c1_cpx*T3(i))+(c2_cpx*(T3(i)**(-2.)))+(c3_cpx*(T3(i)**(-1./2.)))+&
			(c4_cpx*(T3(i)**2.))))+(0.65*(c0_plag+(c1_plag*T3(i))+(c2_plag*(T3(i)**(-2.)))+&
				(c3_plag*(T3(i)**(-1./2.)))+(c4_plag*(T3(i)**2.)))))
    end do
    do i=zc+1,nz
      intalphaT3(i)=(a0*(T3(i)-T0))+((a1/2.)*((T3(i)**2.)-(T0**2.)))
      p3(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT3(i)))
      cp3(i)=(c0_mantle-(c1_mantle*(T3(i)**(-1./2.)))-(c2_mantle*(T3(i)**(-3.))))
    end do    

    cora=-((T2+T3)*(p2*cp2-p1*cp1))/(p3*cp3+p2*cp2)

end subroutine correction

! ------------------------------------------------------------------------
  
subroutine set_var(Ta,cz,czm,p,cp,alphaP0,rhoP0,prdat0,pdat0,dz)
    
    ! calculates the values of key variables
     
    use global
    implicit none
    real(wp), dimension(0:nz), intent(in)  ::  Ta, alphaP0, rhoP0,prdat0,pdat0
    real(wp), dimension(0:nz), intent(out)  ::  cz,czm, p, cp, dz
    real(wp), dimension(0:nz) :: ki, Di, ka, krad
    real(wp), dimension(0:nz) :: dV, intalphaT
    real(wp), dimension(0:nz) :: Dc
    integer  :: i
    
    ! N.B. changes in pressure profile due to extra water filling in accommodation space created by subsidence are 
    ! ignored here --> initial pressure profile is maintained when calculating pressure-dependent density and conductivity
    ! changes.
    ! This avoids the need to iterate over the subsidence and temperature structure multiple times at each timestep to ensure
    ! that extra hydrostatic pressure predicted by the subsidence and subsidence predicted by conductive model with pressure
    ! and temperature-dependent parameters converge and are self-consistent. This can be implemented with a Newton-Raphson
    ! step but slows the code down significantly and makes little difference to the result.
    do i=0,zc
      intalphaT(i)=(a0c*(Ta(i)-T0))+((a1c/2.)*((Ta(i)**2.)-(T0**2.)))
      p(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT(i)))    
      cp(i)=((0.15*((k0fofa)+(k1fofa*(Ta(i)**(-1./2.)))+(k3fofa*(Ta(i)**(-3.)))))+&
	  	(0.2*(c0_cpx+(c1_cpx*Ta(i))+(c2_cpx*(Ta(i)**(-2.)))+(c3_cpx*(Ta(i)**(-1./2.)))+&
			(c4_cpx*(Ta(i)**2.))))+(0.65*(c0_plag+(c1_plag*Ta(i))+(c2_plag*(Ta(i)**(-2.)))+&
				(c3_plag*(Ta(i)**(-1./2.)))+(c4_plag*(Ta(i)**2.)))))
      Dc(i)=(a_crust+(b_crust*exp(-(Ta(i)-T0)/c_crust))+(d_crust*exp(-(Ta(i)-T0)/e_crust)))*1e-6
      krad(i)=Ar*exp(-(((Ta(i)-Tar)**(2.))/(2*(xa**(2.)))))+Br*exp(-(((Ta(i)-Tbr)**(2.))/(2*(xb**(2.)))))
      ka(i)=(p(i)*Dc(i)*cp(i)*exp(dkP*(prdat0(i)/1.e+9)))+krad(i)
      dV(i)=pdat0(i)/p(i)
    end do
    do i=zc+1,nz
      intalphaT(i)=(a0*(Ta(i)-T0))+((a1/2.)*((Ta(i)**2.)-(T0**2.)))
      p(i)=rhoP0(i)*(1.-(alphaP0(i)*intalphaT(i)))
      cp(i)=(c0_mantle-(c1_mantle*(Ta(i)**(-1./2.)))-(c2_mantle*(Ta(i)**(-3.))))
      Dc(i)=(a_ol+(b_ol*exp(-(Ta(i)-T0)/c_ol))+(d_ol*exp(-(Ta(i)-T0)/e_ol)))*1e-6
      krad(i)=Ar*exp(-(((Ta(i)-Tar)**(2.))/(2*(xa**(2.)))))+Br*exp(-(((Ta(i)-Tbr)**(2.))/(2*(xb**(2.)))))
      ka(i)=(p(i)*Dc(i)*cp(i)*exp(dkP*(prdat0(i)/1.e+9)))+krad(i)
      dV(i)=pdat0(i)/p(i)
    end do

    ki(0:nz-1)=(ka(1:nz)+ka(0:nz-1))
    ki(nz)=0.
    dz(0:nz-1)=((dzinit*(dV(1:nz)+dV(0:nz-1)))/2.)
    dz(nz)=dzinit

    Di(0:nz)=ki(0:nz)*(1./2.)

    !include contraction of dz
    cz(0:nz)=Di(0:nz)/dz(0:nz)
    czm(0:nz)=Di(0:nz)/dz(0:nz)

    cz(nz)=0d0
    czm(nz)=0d0

    return
    end subroutine set_var


  ! ------------------------------------------------------------------------
  subroutine set_eqn(Ta,Tn,cz,czm,p,cp,cor,dtb,Tbase,dz)
      
    !    sets up A matrix and b vector to solve Ax=b
    !    advances T(n) to h(n+dt)
    !    calls tridiagonal solver
    
    use global
    implicit none
    real(wp), dimension(0:nz), intent(in)  ::  Ta, cz,czm,cor,p,cp,dz
    real(wp), dimension(0:nz), intent(out)  ::  Tn
    real(wp), dimension(0:nz) :: ab,al,ad,au,x
    real(wp), dimension(0:nz-1) :: beta,dzi
    real(wp), intent(in) :: dtb,Tbase
    integer ::  i
 
    ! include contraction of dz  
    dzi(0:nz-1)=(dz(1:nz)+dz(0:nz-1))/2.
    ! dzi(0) = avg of 0 and 1 points 
    beta(0:nz-1)=dtb/(2.*p(1:nz)*cp(1:nz)*dzi(0:nz-1))
    
    ab(1:nz-1)=Ta(1:nz-1) + beta(0:nz-2)*(cz(1:nz-1)*Ta(2:nz)-(czm(1:nz-1)+cz(0:nz-2))*Ta(1:nz-1)+czm(0:nz-2)&
    *Ta(0:nz-2))+cor(1:nz-1)
    ab(0)=Ta(0) ! ha(0) + beta(0)*(cx(1)*ha(1)-cxm(1)*ha(0)) + dtb/dz
    ab(nz)=Tbase
    
    al(1:nz-1)=-czm(0:nz-2)*beta(0:nz-2)
    ad(1:nz-1)=1d0+(cz(0:nz-2)+czm(1:nz-1))*beta(0:nz-2)
    au(1:nz-1)=-cz(1:nz-1)*beta(0:nz-2)

    al(0)=0._wp
    ad(0)=1d0 ! +czm(0)*beta(0)
    au(0)=0._wp ! -cz(0)*beta(0)


    al(nz)=0. ! -cz(nz)*beta(nz)
    ad(nz)=1d0 ! +czm(nz)*beta(nz)
    au(nz)=0.

    Tn(:)=0.
    call gtri(al,ad,au,ab,x,0,nz)
    do i=0,nz
    if (x(i) .le. 0) then
        Tn(i)=0.
    else 
        Tn(i)=x(i)
    end if
    end do

    return
  end subroutine set_eqn
 
! ------------------------------------------------------------------------
 
  subroutine gtri(l,d,u,b,x,n1,n2)
      
    ! solves tridiagonal system LnXn-1 + DnXn + UnXn+1=Bn (1<=n<=nz)

    implicit none
    integer, intent(in) :: n1, n2    
    real(wp), dimension(n1:n2), intent(in) :: l,d,u,b
    real(wp), dimension(n1:n2), intent(out) :: x
    real(wp), dimension(n1:n2) :: wb
    real(wp) :: t
    integer :: i

    if (d(1) == 0.) then                 ! if diagonal matrix has only 1 value, stop
       write(1,*) 'singular matrix'
       stop
    end if
    x(n1) = d(n1)
    wb(n1) = b(n1)
    do i = n1+1,n2
       t = l(i) / x(i-1)
       x(i) = d(i) - u(i-1)*t
       wb(i) = b(i) - wb(i-1)*t
    end do
    x(n2) = wb(n2) / x(n2)

    do i = n2-1,n1,-1
       x(i) = ( wb(i) - x(i+1)*u(i) ) / x(i)
    end do

  end subroutine gtri
  
! -------------------------------------------------
  
  function brent(xmin,xin)
      
    ! Brent minimisation algorithm to find (V0/V)T 
    
    use global
    implicit none    
    real(wp) :: brent,xin,xmin
    integer, parameter :: itmax=100
    real(wp), parameter :: CGOLD=.3819660
    real(wp), parameter :: ZEPS=1.0e-10
    real(wp), parameter :: AX=1.
    real(wp), parameter :: BX=1.5
    real(wp), parameter :: CX=2.0
    real(wp), parameter :: tol=1.48e-8
    
!     given a function f, and given a bracketing triplet of abscissas AX, bx, cx (such that bx is
!     between AX and cx, and f(bx) is less than both f(AX) and f(CX)), this routine isolates
!     the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
!     the minimum is returned as xmin, and the minimum function value is returned as brent,
!     the returned function value.
!     Parameters: MAXimum allowed number of iterations; golden ratio; and a small number that
!     protects against trying to achieve fractional accuracy for a minimum that happens to be
!     exactly zero.
    integer :: iter
    real(wp) :: a,b,dd,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(AX,CX)        !a and b must be in ascending order, though the input
    b=max(AX,CX)    !abscissas need not be.
    v=BX          !Initializations...
    w=v
    x=v
    e=0.    !This will be the distance moved on the step before last.
    fx=abs((K0*(3./2.)*(x**(7./3.)-x**(5./3.))*(1.+(((3./4.)*(KT-4.))*(x**(2./3.)-1.))))-xin)
    fv=fx
    fw=fx
    do iter=1,itmax   !Main program loop.
!       write(*,*) 'iteration', iter
      xm=0.5*(a+b)
      tol1=(tol*abs(x))+ZEPS
      tol2=2.*tol1
      if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3    !Test for done here.
      if(abs(e).gt.tol1) then        !Construct a trial parabolic fit.
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.*(q-r)
        if(q.gt.0.) p=-p
        q=abs(q)
        etemp=e
        e=dd
    if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
    !The above conditions determine the acceptability of the parabolic fit. Here it is o.k.:
    dd=p/q !Take the parabolic step.
    u=x+dd
    if(u-a.lt.tol2 .or. b-u.lt.tol2) dd=sign(tol1,xm-x)
    goto 2     !Skip over the golden section step.
      endif
1     if(x.ge.xm) then    !We arrive here for a golden section step, which we take
        e=a-x     !into the larger of the two segments.
      else
        e=b-x
      endif
      dd=CGOLD*e           !Take the golden section step.
2     if(abs(dd).ge.tol1) then     !Arrive here with dd computed either from parabolic fit, or
        u=x+dd              !else from golden section.
      else
        u=x+sign(tol1,dd)
      endif
      fu=abs((K0*(3./2.)*(u**(7./3.)-u**(5./3.))*(1.+(((3./4.)*(KT-4.))*(u**(2./3.)-1.))))-xin) !This is the one function evaluation per iteration,
      if(fu.le.fx) then    !and now we have to decide what to do with our function
        if(u.ge.x) then     !evaluation. Housekeeping follows:
          a=x
        else
          b=x
        endif
        v=w
        fv=fw
        w=x
        fw=fx
        x=u
        fx=fu
      else
        if(u.lt.x) then
          a=u
        else
          b=u
        endif
        if(fu.le.fw .or. w.eq.x) then
      v=w
      fv=fw
      w=u
      fw=fu
    else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
      v=u
      fv=fu
    endif
      endif      !Done with housekeeping. Back for another iteration.
    enddo 
    write(*,*) 'brent exceed maximum iterations'
3   xmin=x    !Arrive here ready to exit with best values.
    brent=fx
    return
    end function brent

! -------------------------------------------------

  function rtbis(pb,cont,pw,rd)
      
    ! root bisection to find seafloor subsidence
    
    use global
    implicit none    
    real(wp) :: rtbis,pb,cont,fmid,f,xmid,dx,pw,rd
    real(wp), parameter :: AX=-0.1
    real(wp), parameter :: BX=10.0
    real(wp), parameter :: xacc=1.48e-8
    integer, parameter :: jmax=100
    integer :: j
    
    ! Using bisection, find the root of a function known to lie between AX and BX.
    ! The root, returned as rtbis, will be refined until its accuracy is +/-xacc


    fmid=((pb/(pb-(1028.+(2.4*(BX+rd)))))*cont)-BX
    f=((pb/(pb-(1028.+(2.4*(AX+rd)))))*cont)-AX
    if (f*fmid.ge.0) then
      write(*,*) 'root must be bracketed for bisection'
    endif
    if (f.lt.0) then
      rtbis=AX
      dx=BX-AX
    else
      rtbis=BX
      dx=AX-BX
    endif
    do j=1,jmax
      dx=dx*0.5
      xmid=rtbis+dx
      fmid=((pb/(pb-(1028.+(2.4*(xmid+rd)))))*cont)-xmid
      if (fmid.le.0) then
    rtbis=xmid
      endif
      if (abs(dx).lt.xacc.or.fmid.eq.0) then
          pw=1028.+(2.4*(rtbis+rd))
    return
      endif
    enddo
    write(*,*) 'too many root bisections'
  end function rtbis
    
!--------------------------------------------------

end program RHCW18
