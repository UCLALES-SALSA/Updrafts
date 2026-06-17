!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulator parameters from data files,
! determine emulator input parameters based on column properties,
! and predict updraft velocities.
!
!   Tomi Raatikainen (FMI) 20.3.2026
!
! Current emulator (J25): emulator development as in Nordling et al. (2024) 
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
    CHARACTER(LEN=100) :: fname = 'gp_wp_N24.emu'
    ! 2) Emulator input vector
    REAL, SAVE, PUBLIC :: input_vec(7)
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
        REAL, INTENT(IN) :: input_vec(*) ! Inputs (will be normalized)
        REAL, INTENT(OUT) :: pred ! Predicted unnormalized emulator output
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
        INTEGER :: emu_flag ! Output flag: 0=OK, 1=land, 2=sea ice, 3=fog, 4=no low cloud,
                            ! 5=too high cloud top, 6=ice in the low cloud, 7=clouds above
        !
        ! Run test
        CALL test_column_ifs(klev, pfull, tpot, lwc, iwc, cc, qv, cdnc, lwp, iwp, &
            psrf, lsm, cos_mu, emu_flag, input_vec)
        !
        ! Emulator inputs already set if emu_flag=0
        !
        ! Return flag
        emu_get_inputs = emu_flag
    END FUNCTION emu_get_inputs


    ! Predict updraft velocity for a single column - with some IFS modifications
    SUBROUTINE test_column_ifs(klev, pfull, tpot, lwc, iwc, cc, qv, cdnc, lwp, iwp, psrf, lsm, cos_mu, emu_flag, input)
        INTEGER, INTENT(IN) :: &
            klev            ! Number of vertical levels
        REAL, INTENT(in) :: &
            pfull(klev),  & ! Full-level pressure (Pa)
            tpot(klev),   & ! Potential temperature (K)
            lwc(klev),    & ! Cloud liquid water content (kg/kg)
            iwc(klev),    & ! Cloud ice water content (kg/kg)
            cc(klev),     & ! Cloud cover (-)
            qv(klev),     & ! Specific humidity (kg/kg)
            cdnc(klev),   &  ! CDNC (1/cm3)
            lwp, iwp,     & ! Liquid and ice water paths (kg/m2)
            psrf,         & ! Surface pressure (Pa)
            lsm,          & ! Land cover fraction (-)
            cos_mu          ! Cosine of solar zenith angle (-)
        INTEGER , INTENT(OUT) :: &
            emu_flag        ! Output flag: 0=OK, 1=land, 2=sea ice, 3=fog, 4=no low cloud,
                            ! 5=too high cloud top, 6=ice in the low cloud, 7=clouds above
        REAL, INTENT(OUT) ::  &
            input(7)        ! Emulator input vector 
        ! Local variables
        INTEGER :: i, ind_hi, ind_lo, itop, ibase
        REAL :: low_lwp, low_iwp, p1, p2, dgz, qt(klev)
        REAL, PARAMETER :: &
            g_rcp = 1./9.81, & ! 1/g
            min_lwc = 1e-5, &  ! Cloud water treshold
            min_cc = 0.05,  &  ! Cloud coverage treshold
            pmin = 225e2       ! Maximum cloud top height (p_srfc-p_top)

        emu_flag = 0  ! 0=ok

        ! No land (lsm>0.1) or sea ice (theta(klev)<265.0 K)
        IF (lsm>0.1) THEN
            emu_flag = 1
            RETURN
        ELSEIF (tpot(klev)<265.0) THEN
            emu_flag = 2
            RETURN
        ENDIF

        ! No columns with fog (cloud at the first level above surface)
        !   -Based on minumum cloud water content per cloud area, lwc/cc>min_lwc
        !   -Additional minimum cloud coverage limit, cc>min_cc
        IF (lwc(klev)>min_lwc*cc(klev) .AND. cc(klev)>min_cc) THEN
            emu_flag = 3
            RETURN
        ENDIF

        ! Locate the low cloud base and top levels
        !   -Based on minumum cloud water content per cloud area, lwc/cc>min_lwc
        !   -Additional minimum cloud coverage limit, cc>min_cc
        itop=0
        ibase=0
        ! Calculate LWP and IWP within the low cloud (from surface up to top)
        !   LWP=SUM((rhoa*lwc)*dz)=SUM(lwc*dp/g), where dp=rhoa*g*dz
        low_lwp = 0.0
        low_iwp = 0.0
        p1=psrf ! Surface
        DO i=klev,2,-1
            IF (lwc(i)>min_lwc*cc(i) .AND. cc(i)>min_cc) THEN
                IF (ibase==0) ibase=i
                itop=i
            ELSEIF (itop>0) THEN
                exit
            ENDIF
            !
            p2=0.5*(pfull(i-1)+pfull(i)) ! Interface above
            dgz=(p1-p2)*g_rcp
            IF (cc(i)>min_cc) THEN
                low_lwp=low_lwp+lwc(i)/cc(i)*dgz
                low_iwp=low_iwp+iwc(i)/cc(i)*dgz
            ELSE
                low_lwp=low_lwp+lwc(i)*dgz
                low_iwp=low_iwp+iwc(i)*dgz
            ENDIF
            p1=p2
        ENDDO
        !
        IF (itop==0) THEN
            ! No (low) cloud
            emu_flag = 4
            RETURN
        ELSEIF (psrf-0.5*(pfull(itop-1)+pfull(itop))>pmin) THEN
            ! Too high cloud top height
            emu_flag = 5
            RETURN
        ELSEIF (low_iwp>0.25*low_lwp) THEN
            ! No mixed-phase clouds
            emu_flag = 6
            RETURN
        ELSEIF (low_lwp < (lwp+iwp)*0.5) THEN
            ! No clouds above the low cloud (LWP for the low cloud is at least half of the total LWP+IWP)
            emu_flag = 7
            RETURN
        ENDIF


        ! ******* Emulator inputs for the valid column *******

        ind_lo = MIN(ibase+2,klev)
        ind_hi = itop-2
        ! 1) q_inv (g/kg): inversion strength of total water (vapor and cloud) mass mixing ratio
        ! max-min values of the total water for levels from two layers below cloud up to two levels above cloud
        qt(:) = qv(:) ! Total water for the cloudy part of the column
        DO i=ind_hi,ind_lo
            IF (cc(i)>min_cc) qt(i) = qv(i)+lwc(i)/cc(i)
        ENDDO
        input(1) = 1e3*(MAXVAL(qt(ind_hi:ind_lo))-MINVAL(qt(ind_hi:ind_lo)))
        !
        ! 2) tpot_inv (K): inversion strength of liquid water potential temperature
        ! max-min values of potential temperature (the same levels as for q_inv)
        input(2) = MAXVAL(tpot(ind_hi:ind_lo))-MINVAL(tpot(ind_hi:ind_lo))
        !
        ! 3) lwp (g/m^2): liquid water path for the low cloud
        ! integrated from the surface up to the cloud top
        input(3) = low_lwp*1e3 ! LWP for the cloudy part of the column
        !
        ! 4) tpot_pbl (K): liquid water potential temperature in the boundary layer
        ! min value of potential temperature (the same levels as for q_inv)
        input(4) = MINVAL(tpot(ind_hi:ind_lo))
        !
        ! 5) blh (hPa): inversion/cloud top/planetary boundary layer height described as pressure
        ! difference between surface and cloud top (use half-level/interface pressures)
        input(5) = 0.01*(psrf-0.5*(pfull(itop-1)+pfull(itop)))
        !
        ! 6) CDNC: Cloud droplet number concentration (cdnc) for Level 3 microphysics
        input(6) = SUM(cdnc(itop:ibase))/FLOAT(ibase-itop+1)
        !
        ! 7) cos_mu: cosine of the solar zenith angle
        input(7) = cos_mu

        ! ******* Post-filtering *******
        ! Not applied

    END SUBROUTINE test_column_ifs

END MODULE me_emulator
