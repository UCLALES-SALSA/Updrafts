!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulator parameters from data files,
! determine emulator input parameters based on column properties,
! and predict updraft velocities.
!
!   Tomi Raatikainen (FMI) 20.3.2026
!
! Current emulator (J26): emulator developments from January 2026.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE me_emulator
    USE me_gp, ONLY : DenseGP, predict
    IMPLICIT NONE
    !
    PRIVATE
    PUBLIC :: init_emulators, emu_get_inputs, emu_predict
    !
    ! Emulator
    TYPE(DenseGP), ALLOCATABLE :: gp
    !
    ! Emulator-dependent settings
    ! 1) Default emulator input file
    CHARACTER(LEN=100) :: fname = 'gp_ws_Sc_J26.emu'
    ! 2) Emulator input vector
    REAL, SAVE, PUBLIC :: input_vec(9)

    ! Specific settings

    ! Emulator parameters for cumulus and stratocumulus clouds
    REAL :: lwp_rad, cos_mu_rad ! Common emulator parameters (private)
    REAL :: par_cu(8), par_sc(9) ! Emulator parameters
    REAL :: err_cu, err_sc ! Errors for the linearized temperature and moisture profiles
    INTEGER :: flag_cu, flag_sc ! Flags
    !
    ! Constants
    REAL, PARAMETER :: p00=1.0e5,    & ! Reference pressure (Pa)
                       R=287.04,     & ! Specific gas constant for dry air (R_specific=R/M), J/kg/K
                       Rm=461.5,     & ! -||- for water
                       ep2=Rm/R-1.0, & ! M_air/M_water-1
                       cp=1005.0,    & ! Specific heat for a constant pressure
                       alvl = 2.5e6, & ! Latent heat of vaporization
                       rcp=R/cp 	   ! R/cp
    !
    !
    ! Additional outputs from test_column_ifs and other subroutines for testing and emulator development
    INTEGER :: itop, ibase, imax ! Indexes for cloud top, base and maximum 
    REAL :: low_lwp, low_iwp ! LWP and IWP for the low cloud
    REAL :: err_th_cu, err_rt_cu, par_cu_fit(11) ! Detailed temperature and humidity fits and their errors
    REAL :: err_th_sc, err_rt_sc, par_sc_fit(11)

CONTAINS

    ! Initialize emulators
    SUBROUTINE init_emulators(gp_fname)
        ! Load emulator data file gp_fname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: gp_fname
        !
        ! Use the specified input file
        IF (PRESENT(gp_fname)) fname = gp_fname
        !
        ! Allocate emulator and read the emulator data file
        allocate(gp, source = DenseGP(fname))
        !
    END SUBROUTINE init_emulators

    ! Predict if the inputs are within the range of training data
    SUBROUTINE emu_predict(input_vec,pred,fail)
        REAL, INTENT(INOUT) :: input_vec(*) ! Inputs (will be normalized)
        REAL, INTENT(INOUT) :: pred ! Predicted unnormalized emulator output
        LOGICAL, INTENT(OUT) :: fail ! Flag (OK=.FALSE.)
        INTEGER :: n
        !
        n=SIZE(gp%min_x)
        fail = ANY(input_vec(1:n)<gp%min_x .OR. input_vec(1:n)>gp%max_x)
        IF (fail) THEN
            ! Not within the range of training data
            pred = -999.0
        ELSE
            pred = predict(gp,input_vec(1:n))
        ENDIF
        !
    END SUBROUTINE emu_predict

    ! Calculate emulator inputs for one IFS column
    INTEGER FUNCTION emu_get_inputs(klev, pfull, tpot, lwc, rwc, iwc, cc, qv, cdnc, lwp, iwp, &
            psrf, lhf, shf, lsm, cos_mu)
        INTEGER, INTENT(IN) :: &
            klev            ! Number of vertical levels
        REAL, INTENT(in) :: &
            pfull(klev),  & ! Full-level pressure (Pa)
            tpot(klev),   & ! Potential temperature (K)
            lwc(klev),    & ! Cloud liquid water content (kg/kg)
            rwc(klev),    & ! Rain  water content (kg/kg)
            iwc(klev),    & ! Cloud ice water content (kg/kg)
            cc(klev),     & ! Cloud cover (-)
            qv(klev),     & ! Specific humidity (kg/kg)
            cdnc(klev),   & ! CDNC (1/cm3)
            lwp, iwp,     & ! Liquid and ice water paths (kg/m2)
            psrf,         & ! Surface pressure (Pa)
            lhf, shf,     & ! Latent and sensible heat fluxes (W/m2)
            lsm,          & ! Land cover fraction (-)
            cos_mu          ! Cosine of the solar zenith angle
        ! Variables for the current emulator
        INTEGER :: emu_flag
        !
        ! Run test
        emu_flag = test_column_ifs(klev, pfull, tpot, lwc, iwc, cc, qv, lwp, iwp, psrf, lsm, cos_mu)
        !
        ! Emulator inputs - Sc only
        input_vec(:)=0.0 
        IF (emu_flag==0 .AND. flag_sc==0) THEN
            input_vec(:)=par_sc(:)
            emu_get_inputs=0
        ELSEIF (emu_flag<8) THEN
            emu_get_inputs=emu_flag
        ELSE
            emu_get_inputs=flag_sc
        ENDIF
        !
        ! Return flag
        !emu_get_inputs = flag_sc
    END FUNCTION emu_get_inputs

    ! Determine emulator parameters for an IFS column. Function returns 0 when emulator
    ! parameters have been extracted while other values indicate reason for rejection
    INTEGER FUNCTION test_column_ifs(klev, pfull, tpot, lwc, iwc, cc, qv, &
                        lwp, iwp, psurf, lsm, cos_mu)
        INTEGER, INTENT(IN) :: &
            klev            ! Number of vertical levels
        REAL, INTENT(in) :: &
            pfull(klev),  & ! Full-level pressure (Pa)
            tpot(klev),   & ! Potential temperature (K)
            lwc(klev),    & ! Cloud liquid water content (kg/kg)
            iwc(klev),    & ! Cloud ice water content (kg/kg)
            cc(klev),     & ! Cloud cover (-)
            qv(klev),     & ! Specific humidity (kg/kg)
            lwp, iwp,     & ! Liquid and ice water paths (kg/m2)
            psurf,        & ! Surface pressure (Pa)
            lsm,          & ! Land cover fraction (-)
            cos_mu          ! Cosine of the solar zenith angle
        ! Local variables
        INTEGER :: i, n_clear
        REAL :: p1, p2, dgz
        REAL, PARAMETER :: &
            g_inv = 1./9.81, & ! 1/g
            min_lwc = 1e-5, &  ! Cloud water treshold
            min_cc = 0.05,  &  ! Cloud coverage treshold
            pmin_max = 225e2, & ! Maximum cloud max height (p_surf-p_max)
            pmin_top = 300e2    ! Maximum cloud top height (p_surf-p_top)
        INTEGER, PARAMETER :: &
            n_clear_levels = 2  ! Required number of clear levels above clouds

        ! Remove land and sea ice (theta[srfc]<265 K) columns
        IF (lsm>0.1) THEN
            test_column_ifs = 1
            RETURN
        ELSEIF (tpot(klev)<265.0) THEN
            test_column_ifs = 2
            RETURN
        ENDIF

        ! No columns with fog (cloud at the first level above surface)
        !   -Based on minumum cloud water content per cloud area, lwc/cc>min_lwc
        !   -Additional minimum cloud coverage limit, cc>min_cc
        IF (lwc(klev)>cc(klev)*min_lwc .AND. cc(klev)>min_cc) THEN
            test_column_ifs = 3
            RETURN
        ENDIF

        ! Locate the low cloud base and top levels, and the level with maximum liquid water
        !   -Based on minumum cloud water content per cloud area, lwc/cc>min_lwc
        !   -Additional minimum cloud coverage limit, cc>min_cc
        !   -Find the most common cloud, i.e., the maximum of cc*(lwc/cc)=lwc
        itop = 0
        ibase = 0
        ! Calculate LWP and IWP within the low cloud (from surface up to top)
        !   LWP=SUM((rhoa*lwc)*dz)=SUM(lwc*dp/g), where dp=rhoa*g*dz
        low_lwp = 0.0
        low_iwp = 0.0
        p1 = psurf ! Surface
        n_clear = 0
        DO i=klev,2,-1 ! Start from the last, i.e. surface
            IF (lwc(i)>cc(i)*min_lwc .AND. cc(i)>min_cc) THEN
                IF (ibase==0) THEN
                    ibase=i
                    imax=i ! The maximum cloud water
                ENDIF
                itop=i
                IF (lwc(i)>lwc(imax)) imax=i ! The most common top
                n_clear = 0
            ELSEIF (itop>0) THEN
                n_clear=n_clear+1 ! This is clear
                IF (n_clear==n_clear_levels) exit ! Stop when found enough clear levels above cloud
            ENDIF
            !
            p2=0.5*(pfull(i-1)+pfull(i)) ! Interface above
            dgz=(p1-p2)*g_inv
            low_lwp=low_lwp+lwc(i)*dgz
            low_iwp=low_iwp+iwc(i)*dgz
            p1=p2
        ENDDO
        !
        IF (itop==0) THEN
            ! Low cloud not found
            test_column_ifs = 4
            RETURN
        ELSEIF (psurf-pfull(imax)>pmin_max .OR. psurf-pfull(itop)>pmin_top) THEN
            ! Too high cloud top height
            test_column_ifs = 5
            RETURN
        ELSEIF (low_iwp>0.25*low_lwp) THEN
            ! No mixed-phase clouds
            test_column_ifs = 6
            RETURN
        ENDIF

        ! Obtain and test common variables
        ! 1) LWP_rad (g/m2): LWP and IWP above the cloud for calculating radiative fluxes
        ! - Liquid is has about six times larger impact on optical depth than ice,
        !   so add ice divided by 6 and consider all as liquid
        lwp_rad=(MAX(0.0,(lwp-low_lwp)) + MAX(0.0,(iwp-low_iwp))/6.0)*1e3
        IF (lwp_rad>400.0) THEN
            ! LWP above 400 g/m2
            test_column_ifs = 7
            RETURN
        ENDIF
        !
        ! 2)  Cosine of the solar zenith angle (any value is good)
        cos_mu_rad = cos_mu

        ! Parameterize cloud temperature and humidity profiles (returns zero if successful)
        CALL cloud_par_simple(klev,pfull,tpot,qv,lwc,cc,psurf,ibase,imax,itop)
        IF (flag_sc/=0 .AND. flag_cu/=0) THEN
            ! No solution because both Sc and Cu parameterizations failed
            test_column_ifs = 10*flag_cu+flag_sc
        ELSE
            ! At least one parameterization is good
            test_column_ifs = 0
        ENDIF

        ! All done
        RETURN
    END FUNCTION test_column_ifs

    SUBROUTINE cloud_par_simple(n,press,theta,rv,rc,cc,p_surf,ibase,imax,itop)
        INTEGER, INTENT(IN) :: n   ! Number of vertical levels
        REAL, INTENT(in) :: press(n),theta(n),rv(n),rc(n),cc(n),p_surf
        INTEGER, INTENT(in) :: ibase, imax, itop
        ! Local
        REAL :: rc_cc(n), theta_L(n), rt(n)
        !REAL :: err_th_cu, err_rt_cu, par_cu_fit(11)
        !REAL :: err_th_sc, err_rt_sc, par_sc_fit(11)
        INTEGER :: first
        !
        ! Not all model levels are relevant: include cloud maximum+4 (and cloud top if above that)
        first=min(imax-4,itop) ! imax-4>1
        !
        ! Cumulus
        ! ======
        ! Focus on the whole column, which is not fully cloudy
        err_cu=1e10
        err_rt_cu=1e10
        flag_cu=9 ! If theta fit fails
        CALL cu_par_th(first,n,press,theta,p_surf,ibase,err_th_cu,par_cu_fit)
        IF (err_th_cu<1.0) THEN
            rt(:) = rv(:) + rc(:) ! Total water
            CALL cu_par_rt(first,n,press,rt,err_rt_cu,par_cu_fit)
            IF (err_rt_cu<1.0) THEN
                ! Simplified cumulus parameters (flag=0 means valid results)
                CALL simplify_cu_pars(par_cu_fit,flag_cu)
                if (flag_cu==0) err_cu=SQRT(err_th_cu**2+err_rt_cu**2)
            ELSE
                flag_cu=8
            ENDIF
        ENDIF
        !
        !
        ! Stratocumulus
        ! ==========
        ! Focus on the cloudy portion of the column
        !
        ! Derived thermodynamic variables for the cloudy portion of the column
        ! a) LWC
        rc_cc=rc
        WHERE (cc>0.05) rc_cc=rc/cc
        ! b) Liquid water potential temperature
        theta_L(:) = theta(:)-rc_cc(:)*(alvl/cp)/(press(:)/p00)**rcp
        !
        err_sc=1e10
        err_rt_sc=1e10
        flag_sc=9
        CALL sc_par_th(first,n,press,theta_L,p_surf,ibase,imax,err_th_sc,par_sc_fit)
        IF (err_th_sc<1.0) THEN
            rt(:) = rv(:) + rc_cc(:) ! Total water
            CALL sc_par_rt(first,n,press,rt,err_rt_sc,par_sc_fit)
            IF (err_rt_sc<1.0) THEN
                ! Simplified stratocumulus parameters
                CALL simplify_sc_pars(par_sc_fit,rc_cc(imax),press(imax),p_surf,flag_sc)
                if (flag_sc==0) err_sc=SQRT(err_th_sc**2+err_rt_sc**2)
            ELSE
                flag_sc=8
            ENDIF
        ENDIF
        !
    END SUBROUTINE cloud_par_simple


    ! ************** Cumulus parameterization **************

    SUBROUTINE cu_par_th(first,n,press,theta,p_surf,ibase,err_th,par)
        ! Describe cumulus cloud cloud potential temperature profiles with linear sections
        INTEGER, INTENT(IN) :: first,n ! The first and last vertical level
        REAL, INTENT(in) :: press(n),theta(n),p_surf
        INTEGER, INTENT(in) :: ibase
        REAL, INTENT(out) :: err_th,par(11)
        ! Local parameters
        INTEGER :: i, m
        REAL :: th_srfc,psrfc,th_mix,pmix,th_top,ptop
        REAL :: sx, sy, sxy, sx2
        REAL :: a, b, err
        
        REAL :: tavg, q_sat, lapse
        !
        err_th=1e10
        par(:)=0.0
        !
        m=0
        sx=0.0; sy=0.0; sxy=0.0; sx2=0.0
        DO i=first,n-1
            ! Linear section from top to point i
            sx=sx+press(i)
            sy=sy+theta(i)
            sxy=sxy+press(i)*theta(i)
            sx2=sx2+press(i)**2
            m=m+1
            !
            ! Linear inversions
            if (i>=ibase) THEN
                ! Linear (cumulus inversion) from first to i, which includes the cloud
                ! Well-mixed boundary layer from i+1 to n
                !
                ! Linear section from first to i
                b=(m*sxy-sx*sy)/(m*sx2-sx**2)
                a=sy/m-b*sx/m
                !
                ! dtheta/dp must be negative (about -0.0005 K/Pa)
                !if (b>-1e-4) CYCLE
                !
                ! Single well-mixed boundary layer from i+1 to n
                th_srfc=SUM(theta(i+1:n))/(n-i)
                psrfc=p_surf
                !
                ! Additional constraint for the lapse rate (b) based on moist adiabatic lapse rate
                tavg=th_srfc*(press(i)/p00)**rcp ! Temperature
                q_sat=rslf(press(i),tavg) ! Saturation mixing ratio
                ! Lapse rate, d(theta)/dp
                lapse=-R*th_srfc/press(i)*( 1/cp - (1+alvl*q_sat/(R*tavg))/(cp+alvl**2*q_sat/(Rm*tavg**2)) )
                ! Conditionally unstable when the lapse rate is larger than the dry adiabatic
                ! lapse rate (d(theta)/dp=0) but smaller than the moist adiabatic lapse rate.
                if (b<lapse*1.5 .OR. b>0.25*lapse) CYCLE ! Here dtheta/dp<0 while dT/dz is positive
                !
                ! Need to have continuous theta, so find the intersect between constant th_srfc and the linear section:
                pmix=(th_srfc-a)/b	! th_srfc=a+b*pmix => pmix=...
                ! Adjust boundary layer theta if the intersect pmix is not press(i)<=pmix<press(i+1)
                if (pmix<press(i)) then
                    ! Boundary layer theta is higher than the value based on the linear section
                    pmix=press(i)
                    th_srfc=a+b*pmix
                elseif (pmix>press(i+1)*0.999) then
                    ! Boundary layer theta is lower.
                    pmix=press(i+1)*0.999
                    th_srfc=a+b*pmix
                endif
                !
                ! Here top = mixed; slope b is used above mixed layer
                ptop=pmix; th_mix=th_srfc; th_top=th_srfc
                !
                !  Error calculations
                err=calc_error(first,n,press,theta,psrfc,th_srfc,pmix,th_mix,ptop,th_top,b)
                if (err<err_th) then
                    par=(/psrfc, th_srfc, pmix, th_mix, ptop, th_top, 0.,0.,0.0, b,0.0/)
                    err_th=err
                endif
            endif
        ENDDO
        !
    end SUBROUTINE cu_par_th
    !
    SUBROUTINE cu_par_rt(first,n,press,rt,err_rt,par)
        ! Describe cumulus cloud cloud total water mixing ratio profiles with linear sections
        INTEGER, INTENT(IN) :: first,n   ! The first and last vertical level
        REAL, INTENT(in) :: press(n),rt(n)
        REAL, INTENT(out) :: err_rt
        REAL, INTENT(inout) :: par(11)
        ! Local
        INTEGER :: i, m, n_mix
        REAL :: rt_srfc,psrfc,rt_mix,pmix,rt_top,ptop,drtdp
        REAL :: sx, sy, sxy, sx2, a, b
        !
        ! Pressure levels
        psrfc=par(1)
        pmix=par(3)
        ptop=par(5) ! ptop=pmix
        !
        ! Parameterize total water
        m=0
        sx=0.0; sy=0.0; sxy=0.0; sx2=0.0
        rt_mix=0.0; n_mix=0
        DO i=first,n
            IF (press(i)<=ptop) THEN
                ! Linear section from top to ptop
                sx=sx+press(i)
                sy=sy+rt(i)
                sxy=sxy+press(i)*rt(i)
                sx2=sx2+press(i)**2
                m=m+1
            ELSE
                ! Constant boundary layer
   				rt_mix=rt_mix+rt(i)
				n_mix=n_mix+1
            ENDIF
        ENDDO
        !
		! Linear above top (m>0)
		b=max(0.0, (m*sxy-sx*sy)/(m*sx2-sx**2) ) ! Must be positive
		a=sy/m-b*sx/m
		!
		rt_top=(a+b*ptop) ! kg/kg
		drtdp=b 		  ! kg/kg/Pa
        !
        ! Mixed layer (n_mix>0): the linear section is important for cumulus clouds,
        ! so allow humidity jump below the linear section
        rt_mix=rt_mix/n_mix
        !
        ! No de-coupling
        rt_srfc=rt_mix
        !
        ! Save/update
		par(7)=rt_top
		par(8)=rt_mix
		par(9)=rt_srfc
		par(11)=drtdp
         !
        ! Calculate error (1e6 is for conversion from (kg/kg)^2 to (g/kg)^2)
        err_rt= calc_error(first,n,press,rt,psrfc,rt_srfc,pmix,rt_mix,ptop,rt_top,drtdp)*1e6
        !
    end SUBROUTINE cu_par_rt
    !
    SUBROUTINE simplify_cu_pars(par,flag)
        ! Simplify cumulus cloud temperature and humidity profile descriptions for emulators
        REAL, INTENT(in) :: par(11)
        INTEGER, INTENT(OUT) :: flag
        REAL :: tavg, rs
        !
        !             1          2              3        4              5        6            7            8            9            10         11
        ! par=(/psrfc, th_srfc, pmix, th_mix, ptop, th_top, rt_top, rt_mix, rt_srfc, dthdp, drtdp/)
        ! Here top = mixed, so ptop=pmix and th_top=th_mix=th_srfc
        !
        ! Surface layer (liquid water) potential temperature (K)
        par_cu(1) = par(2)
        ! Saturation ratio at the top of surface layer (-)
        tavg=par(2)*(par(5)/p00)**rcp   ! Temperature (K)
        rs=rslf(par(5),tavg)            ! Saturation mixing ratio (kg/kg)  
        par_cu(2) = min(par(7),par(9))/rs
        !
        ! Total water mixing ratio inversion (kg/kg)
        par_cu(3) = par(7)-par(9)
        ! Boundary layer thickness (Pa)
        par_cu(4) = par(1)-par(5)
        !
        ! Slopes
        par_cu(5) = par(10) ! dth/dp (K/Pa)
        par_cu(6) = par(11) ! drt/dp (kg/kg/Pa)
        !
        ! Common parameters
        par_cu(7) = lwp_rad
        par_cu(8) = cos_mu_rad
        !
        ! Test if these are reasonble
        flag=0
        IF (par_cu(1)<265.0 .OR. par_cu(1)>305.0) THEN
            ! Surface layer (liquid water) potential temperature (K)
            flag=1
        ELSEIF (par_cu(2)<0.5 .OR. par_cu(2)>1.1) THEN
            ! Saturation ratio at the top of surface layer (-)
            flag=2
        ELSEIF (par_cu(3)<-4e-3 .OR. par_cu(3)>2e-3) THEN
            ! Total water mixing ratio inversion (kg/kg)
            flag=3
        ELSEIF (par_cu(4)<200.0 .OR. par_cu(4)>15e3) THEN
            ! Boundary layer thickness (Pa)
            flag=4
        ELSEIF (par_cu(5)<-2e-3 .OR. par_cu(5)>-1e-4) THEN
            ! dth/dp (K/Pa)
            flag=5
        ELSEIF (par_cu(6)<1e-8 .OR. par_cu(6)>1.5e-6) THEN
            !  drt/dp (kg/kg/Pa)
            flag=6
        ENDIF
        !
    END SUBROUTINE simplify_cu_pars


    ! ************** Stratocumulus parameterization **************

    SUBROUTINE sc_par_th(first,n,press,theta,p_surf,ibase,imax,err_th,par)
        ! Describe stratocumulus cloud potential temperature profiles with linear sections
        INTEGER, INTENT(IN) :: first,n   ! The first and last vertical level
        REAL, INTENT(in) :: press(n),theta(n),p_surf
        INTEGER, INTENT(in) :: ibase, imax
        REAL, INTENT(out) :: err_th, par(11)
        ! Local
        INTEGER :: i, m
        REAL :: th_srfc,psrfc,th_mix,pmix,th_top,ptop
        REAL :: sx, sy, sxy, sx2
        REAL :: a, b, a2, b2, err
        !
        err_th=1e10
        par(:)=0.0
        !
        ! Statistics for the default boundary layer from imax to surface: th_srfc, psrfc, th_mix
        call decoupled_boundary_layer(n,press,theta,p_surf,imax,ibase,th_srfc,psrfc,th_mix)
        !
        sx=0.0; sy=0.0; sxy=0.0; sx2=0.0; m=0
        DO i=first,imax-1
            ! Linear section from top to point i
            sx=sx+press(i)
            sy=sy+theta(i)
            sxy=sxy+press(i)*theta(i)
            sx2=sx2+press(i)**2
            m=m+1
            !
            ! Steep inversions above cloud
            if (i==imax-2) then
                ! Linear (free troposphere) from 1 to i=imax-2
                ! Steep inversion jump at single point just above the cloud (imax-1)
                ! Possibly decoupled boundary layer from imax to n (surface) includes
                ! the cloud; parameters: th_srfc, psrfc, th_mix
                !
                ! Linear section from 0 to i=imax-2
                b=min(0.0,(m*sxy-sx*sy)/(m*sx2-sx**2)) ! Positive values not allowed
                a=(sy-b*sx)/m
                !
                ! Steep inversion jump at single point imax-1 (just above the cloud)
                !   theta(imax-2)=a+b*press(imax-2), theta(imax-1)=theta(imax-1), theta(imax)=th_mix
                !   Take the steeper slope and increase it by 10 %
                b2=1.1*min( (a+b*press(imax-2)-theta(imax-1))/(press(imax-2)-press(imax-1)), &
                    (th_mix-theta(imax-1))/(press(imax)-press(imax-1)) )
                a2=theta(imax-1)-b2*press(imax-1) ! a: y=(y0-k*x0)+k*x
                !
                ! The top of inversion: find the pressure where the two linear profiles intersect
                ! 	y above: y=a+b*x
                ! 	y below: y=a2+b2*x
                ptop=(a2-a)/(b-b2)
                ! Typically this is close to press(imax-1). A slightly larger value would mean that
                ! level imax-1 is in free troposphere, which will be tested next, A value less than
                ! press[imax-2] would mean that level imax-1 is in boundary layer, but such
                ! a large difference is unexpected
                !
                ! Calculate temperature based on the linear free troposphere.
                th_top=a+b*ptop
                !
                if (press(imax-1)>ptop .AND. ptop>press(imax-2) .AND. th_mix<=th_top) then
                    ! Need to have continuous theta between the boundary layer and linear inversion, so find
                    ! the intersect between constant th_mix and the steep inversion:
                    pmix=(th_mix-a2)/b2 ! th_mix=a+b*pmix => pmix=...
                    ! Typically this is close to press[imax-1].
                    !
                    ! Error calculations
                    err=calc_error(first,n,press,theta,psrfc,th_srfc,pmix,th_mix,ptop,th_top,b)
                    if (err<err_th) then
                        par=(/psrfc, th_srfc, pmix, th_mix, ptop, th_top, 0.0,0.0,0.0, b,0.0/)
                        err_th=err
                    endif
                endif
             elseif (i==imax-1) then
                ! Linear (free tropospehere) from 1 to i=imax-1
                ! Steep inversion jump takes place half way between points imax and imax-1
                ! Possibly decoupled boundary layer from imax to n (surface) includes the cloud;
                ! parameters: th_srfc, psrfc, th_mix
                !
                ! Linear section from 0 to i=imax-1
                b=min(0.0,(m*sxy-sx*sy)/(m*sx2-sx**2)) ! Positive values not allowed
                a=(sy-b*sx)/m
                !
                ! Need to have continuous theta, so find the intersect between constant th_mix and the linear section:
                if (b==0.0 .OR. sy/m<th_mix) then
                    ! No intersect (b=0) or the average above is less than that below
                    th_top=max(a,th_mix)
                    ptop=0.01*press(imax)+0.99*press(imax-1) ! Closer to the top
                    pmix=0.99*press(imax)+0.01*press(imax-1) ! Closer to the cloud
                else
                    pmix=(th_mix-a)/b ! th_mix=a+b*pmix => pmix=...
                    if (press(imax-1)<pmix .AND. pmix<press(imax)) then
                        ! The intercept is between press[imax-1] and press[imax], so no inversion
                        ptop=pmix
                        th_top=th_mix !a+b*ptop
                    else
                        ! Inversion between between imax-1 and imax
                        ptop=0.01*press(imax)+0.99*press(imax-1) ! Closer to the top
                        th_top=a+b*ptop ! Based on free tropospehere
                        pmix=0.99*press(imax)+0.01*press(imax-1) ! Closer to the cloud
                    endif
                endif
                !
                ! Boundary layer temperature must be lower than that of free troposphere
                if (th_mix<=th_top) then
                    ! Error calculations
                    err=calc_error(first,n,press,theta,psrfc,th_srfc,pmix,th_mix,ptop,th_top,b)
                    if (err<err_th) then
                        par=(/psrfc, th_srfc, pmix, th_mix, ptop, th_top, 0.0,0.0,0.0, b,0.0/)
                        err_th=err
                    endif
                endif
            endif
        enddo
        !
    END SUBROUTINE sc_par_th
    !
    SUBROUTINE sc_par_rt(first,n,press,rt,err_rt,par)
        ! Describe stratocumulus cloud total water mixing ratio profiles with linear sections
        INTEGER, INTENT(IN) :: first,n   ! The first and last vertical level
        REAL, INTENT(in) :: press(n),rt(n)
        REAL, INTENT(out) :: err_rt
        REAL, INTENT(inout) :: par(11)
        ! Local
        INTEGER :: i, m, n_mix, n_srfc
        REAL :: rt_srfc,psrfc,rt_mix,pmix,rt_top,ptop,drtdp
        REAL :: sx, sy, sxy, sx2, a, b
        !
        ! Pressure levels
        psrfc=par(1)
        pmix=par(3)
        ptop=par(5)
        !
        ! Parameterize total water
		rt_mix=0.0; rt_srfc=0.0
		n_mix=0; n_srfc=0
		sx=0.0; sy=0.0; sxy=0.0; sx2=0.0; m=0
		DO i=first,n
			if (press(i)<=ptop) then
				! Free troposphere: linear
				sx=sx+press(i)
				sy=sy+rt(i)
				sxy=sxy+press(i)*rt(i)
				sx2=sx2+press(i)**2
				m=m+1
    		elseif (press(i)<pmix) then
				! Sc inversion layer: a jump between free troposphere and boundary layer
				!mi=mi+1
        	elseif (press(i)<psrfc) then
				! Cloud layer: fixed
				rt_mix=rt_mix+rt(i)
				n_mix=n_mix+1
            else
				! Surface layer: fixed
				rt_srfc=rt_srfc+rt(i)
				n_srfc=n_srfc+1
			endif
		enddo
		!
        ! Free troposphere: linear
		b=max(0.0, (m*sxy-sx*sy)/(m*sx2-sx**2) ) ! Must be positive
		a=sy/m-b*sx/m
		!
		rt_top=(a+b*ptop) ! kg/kg
		drtdp=b 			! kg/kg/Pa
		!
		! Mixed layer
		rt_mix=rt_mix/n_mix
		!
		! Top and mixed layer need to be continuous and reasonble
		! Boundary layer cloud, so rt_mix represents the cloud => adjust the linear section
		if (ptop==pmix .OR. rt_top>rt_mix) then ! No inversion or the mixed layer is more humid
			! A line that goes through point (pmix,rt_mix):
			drtdp=max(0.0, (sy-m*rt_mix)/(sx-m*pmix) ) ! y(p)=rt_mix+drtdp*(p-pmix)
			! Calculate rt_top=y(ptop)=...
			rt_top=rt_mix+drtdp*(ptop-pmix)
		endif
		!
    	! Surface
		if (n_srfc>0) then
			rt_srfc=rt_srfc/n_srfc
		else
			! No surface layer
			rt_srfc=rt_mix
		endif
		!
		! Save/update
		par(7)=rt_top
		par(8)=rt_mix
		par(9)=rt_srfc
		par(11)=drtdp
        !
        ! Calculate error (1e6 is for conversion from (kg/kg)^2 to (g/kg)^2)
        err_rt = calc_error(first,n,press,rt,psrfc,rt_srfc,pmix,rt_mix,ptop,rt_top,drtdp)*1e6
        !
    end SUBROUTINE sc_par_rt
    !
    SUBROUTINE simplify_sc_pars(par,clw,press,p_surf,flag)
        ! Simplify stratocumulus cloud temperature and humidity profile descriptions for emulators
        REAL, INTENT(in) :: par(11), clw, press, p_surf
        INTEGER, INTENT(OUT) :: flag
        REAL :: dp, tavg, rs, rt_mix, xsi, clw_max, s
        !
        !             1          2             3         4             5        6             7            8            9            10          11
        ! par=(/psrfc, th_srfc, pmix, th_mix, ptop, th_top, rt_top, rt_mix, rt_srfc, dthdp, drtdp/)
        !
        ! Include at least 2.5 kPa (about one model level) for the jumps
        dp=max(0.0, 2.5e3-(par(3)-par(5)) )
        !
        ! Calculate cloud layer total water mixing ratio from the given cloud water: total=rs(press)+clw(press)
        tavg=par(4)*(press/p00)**rcp+clw*alvl/cp ! Temperature (K)
        rs=rslf(press,tavg) ! Saturation mixing ratio (-)
        rt_mix=rs+clw ! Total water instead of par(8)
        !
        ! Fitted cloud top, par(3), may not be the same as "press", so calculate cloud water at the actual cloud top
        ! Moist adiabatic lapse rate
        xsi=(1+ep2*rs)
        tavg=tavg-(1+alvl*rs/(R*tavg))/(cp+alvl**2*rs/(Rm*tavg**2))*(R*tavg*xsi)/press*(press-par(3))
        ! Saturation mixing ratio at the cloud top
        rs=rslf(par(3),tavg)
        ! Cloud water
        clw_max=rt_mix-rs
        !
        ! ** Cloudy layer **
        ! Cloud layer (liquid water) potential temperature (K)
        par_sc(1) = par(4)
        ! Maximum cloud water mixing ratio (kg/kg)
        par_sc(2) = clw_max
        ! Boundary layer thickness (Pa)
        par_sc(3) = p_surf-par(3)
        ! Inversions (dthdp<0 and drtdp>0)
        par_sc(4) = par(6)-par(4)-par(10)*dp ! d(theta)/dp: jump + linear (positive)
		par_sc(5) = par(7)-rt_mix-par(11)*dp ! d(rt)/dp: jump + linear (positive)
        !
        ! ** Decoupling **
        ! Surface layer potential temperature difference from that of the cloud layer (K)
        par_sc(6) = par(4)  - par(2)
        ! Surface layer total water mixing ratio difference from that of the cloud layer (kg/kg)
        par_sc(7) = rt_mix - par(9)
        !
        ! Common parameters
        par_sc(8) = lwp_rad
        par_sc(9) = cos_mu_rad
        !
        ! Calculated and fitted cloud total water mixing ratios should be similar
        s=(rt_mix-par(8))/rt_mix
        !
        ! Test if these are reasonble (cloud only)
        flag=0
        IF (.TRUE.) THEN
            ! Not tested
        ELSEIF (par_sc(1)<265.0 .OR. par_sc(1)>305.0) THEN
            ! Cloud layer (liquid water) potential temperature (K)
            flag=1
        ELSEIF (par_sc(2)<0.05e-3 .OR. par_sc(2)>1.5e-3) THEN
            ! Maximum cloud water mixing ratio (kg/kg)
            flag=2
        ELSEIF (par_sc(3)<2.5e3 .OR. par_sc(3)>22e3) THEN
            ! Boundary layer thickness (Pa)
            flag=3
        ELSEIF (par_sc(4)<1.0+alvl/cp*clw_max .OR. par_sc(4)>20.0) THEN
            ! Potential temperature inversion (K)
            flag=4
        ELSEIF (par_sc(5)<-15e3 .OR. par_sc(5)>-0.05e-3) THEN
            ! Total water mixing ratio inversion (kg/kg)
            flag=5
        ELSEIF (s<-0.15 .OR. s>0.1) THEN
            ! Adjustment to total water
            flag=7
        ENDIF
        !
    END SUBROUTINE simplify_sc_pars
    


    ! Common functions


    ! Prediction error for the linear model
    REAL FUNCTION calc_error(first,n,press,theta,psrfc,ysrfc,pmix,ymix,ptop,ytop,dthdp)
        INTEGER, INTENT(IN) :: first,n   ! Number of vertical levels
        REAL, INTENT(in) :: press(n),theta(n)
        REAL, INTENT(IN) :: psrfc,ysrfc,pmix,ymix,ptop,ytop,dthdp
        INTEGER :: i
        REAL :: err, pred
        !
        ! Error
        err=0.0
        do i=first,n
            ! Prediction
            if (press(i)<=ptop) then
                ! Free troposphere: fixed gradient
                pred=ytop+dthdp*(press(i)-ptop)
            elseif (press(i)<pmix) then
                ! Inversion layer: linear
                pred=ymix+(ytop-ymix)/(ptop-pmix)*(press(i)-pmix)
            elseif (press(i)<psrfc) then
                ! Cloud layer: fixed
                pred=ymix
            else
                ! Surface layer: fixed
                pred=ysrfc
            endif
            ! Error
            err=err+(pred-theta(i))**2
        enddo
        calc_error=err/(n-first+1)
    end FUNCTION calc_error

    ! De-coupled bounday layer for the linear model
    SUBROUTINE decoupled_boundary_layer(n,press,theta,p_surf,itop,ibase,ysrfc,psrfc,ymix)
        INTEGER, INTENT(IN) :: n   ! Number of vertical levels
        REAL, INTENT(in) :: press(n),theta(n),p_surf
        INTEGER, INTENT(in) :: itop, ibase ! Cloudy from top (itop) to base (ibase), which is well mixed
        REAL, INTENT(out) ::ysrfc,psrfc,ymix ! Output parameters descibing the levels with constant theta
        ! Local
        INTEGER :: i, ni, nj
        REAL ::  syi, syi2, syj, syj2, err, err_min
        !
        ! Calculate the mean theta - the default for the whole boundary layer
        syi=0.0; syi2=0.0; ni=n-itop+1
        DO i=itop,n
            syi=syi+theta(i)
            syi2=syi2+theta(i)**2
        enddo
        !
        ! If just one mixed layer
        ymix=syi/ni
        ysrfc=ymix
        psrfc=p_surf
        !
        if (ni>=4) then
            syj=0.0; syj2=0.0; nj=0
            err_min=syi2-syi**2/ni
            !
            ! Split into two regions
            do i=itop,n-1
                ! Points from itop to i (mixed cloud layer)
                syj=syj+theta(i)
                syj2=syj2+theta(i)**2
                nj=nj+1
                ! From i+1 to n (surface)
                syi=syi-theta(i)
                syi2=syi2-theta(i)**2
                ni=ni-1
                ! Sum of squared errors
                err=syi2-syi**2/ni+syj2-syj**2/nj
                ! Require at least two points for each, the cloud is fully in the well-mixed layer, and that surface is cooler
                if (err<err_min .AND. nj>=2 .AND. ni>=2 .AND. i>ibase .AND. syj/nj>syi/ni) then
                    err_min=err
                    ysrfc=syi/ni
                    ymix=syj/nj
                    psrfc=0.5*(press(i)+press(i+1)) ! Half-way
                endif
            enddo
        endif
    end SUBROUTINE decoupled_boundary_layer


! From UCLALES-SALSA thrm.f90
!
! ---------------------------------------------------------------------
! This function calculates the water saturation vapor mixing ratio as
! a function of temperature and pressure
!
  real function rslf(p,t)
  real, intent (in) :: p, t
  real ::  e
  e=esl(t)
  rslf=.622*e/(p-e)
  end function rslf

  real elemental function esl(t)
  real, intent (in) :: t
  real, parameter :: c0=0.6105851e+03, c1=0.4440316e+02,    &
                     c2=0.1430341e+01, c3=0.2641412e-01,    &
                     c4=0.2995057e-03, c5=0.2031998e-05,    &
                     c6=0.6936113e-08, c7=0.2564861e-11,    &
                     c8=-.3704404e-13
  real  :: x
  x=min(max(-80.,t-273.16),50.)
  esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  end function esl


END MODULE me_emulator
