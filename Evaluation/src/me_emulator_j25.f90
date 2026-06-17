!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulator parameters from data files,
! determine emulator input parameters based on column properties,
! and predict updraft velocities.
!
!   Tomi Raatikainen (FMI) 20.3.2026
!
! Current emulator (J25): emulator developments from January 2025.
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
    CHARACTER(LEN=100) :: fname = 'gp_wp_J25.emu'
    ! 2) Emulator input vector
    REAL, SAVE, PUBLIC :: input_vec(9)
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
                            ! 5=too high cloud top, 6=ice in the low cloud, 7=too much ice
        REAL :: q_inv, tpot_inv, tpot_pbl, pbl, lwc_max, cdnc_max, lwp_rad, lhf_out ! Emulator inputs
        !
        ! Run test
        CALL test_column_ifs(klev, pfull, tpot, lwc, rwc, iwc, cc, qv, cdnc, lwp, iwp, &
            psrf, lhf, lsm, emu_flag, q_inv, tpot_inv, tpot_pbl, pbl, lwc_max, cdnc_max, lwp_rad, lhf_out)
        !
        ! Emulator inputs
        input_vec(:)=0.0 
        IF (emu_flag==0) &
            input_vec=(/q_inv, tpot_inv, tpot_pbl, pbl, lwc_max, cdnc_max, cos_mu, lwp_rad, lhf_out/)
        !
        ! Return flag
        emu_get_inputs = emu_flag
    END FUNCTION emu_get_inputs


    ! Calculate emulator variables for one IFS column; originally from process_sample.f90
    SUBROUTINE test_column_ifs(klev, pfull, tpot, lwc, rwc, iwc, cc, qv, cdnc, lwp, iwp, &
            psrf, lhf_ifs, lsm, emu_flag, q_inv, tpot_inv, tpot_pbl, pbl, lwc_max, cdnc_max, lwp_rad, lhf)
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
            lhf_ifs,      & ! Latent heat flux (W/m2)
            lsm             ! Land cover fraction (-)
        INTEGER, INTENT(OUT) :: &
            emu_flag        ! Output flag: 0=OK, 1=land, 2=sea ice, 3=fog, 4=no low cloud,
                            ! 5=too high cloud top, 6=ice in the low cloud, 7=too much ice

        REAL, INTENT(OUT) :: &
            q_inv, tpot_inv, tpot_pbl, pbl, lwc_max, cdnc_max, lwp_rad, lhf ! Emulator inputs
                            
        ! Local variables
        INTEGER :: i, ind_hi, ind_lo, itop, ibase, imax
        REAL :: p1, p2, dgz, low_lwp, low_iwp
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
        IF (lwc(klev)>cc(klev)*min_lwc .AND. cc(klev)>min_cc) THEN
            emu_flag = 3
            RETURN
        ENDIF

        ! Locate the low cloud base and top levels, and the level with maximum liquid water
        !   -Based on minumum cloud water content per cloud area, lwc/cc>min_lwc
        !   -Additional minimum cloud coverage limit, cc>min_cc
        !   -Find the most common cloud, i.e., the maximum of cc*(lwc/cc)=lwc
        itop=0
        ibase=0
        ! Calculate LWP and IWP within the low cloud (from surface up to top)
        !   LWP=SUM((rhoa*lwc)*dz)=SUM(lwc*dp/g), where dp=rhoa*g*dz
        low_lwp = 0.0
        low_iwp = 0.0
        p1=psrf ! Surface
        DO i=klev,2,-1 ! Start from the last, i.e. surface
            IF (lwc(i)>cc(i)*min_lwc .AND. cc(i)>min_cc) THEN
                IF (ibase==0) THEN
                    ibase=i
                    imax=i ! The maximum cloud water
                ENDIF
                itop=i
                IF (lwc(i)>lwc(imax)) imax=i ! The most common top
            ELSEIF (itop>0) THEN
                exit
            ENDIF
            !
            p2=0.5*(pfull(i-1)+pfull(i)) ! Interface above
            dgz=(p1-p2)*g_rcp
            low_lwp=low_lwp+lwc(i)*dgz
            low_iwp=low_iwp+iwc(i)*dgz
            p1=p2
        ENDDO
        !
        IF (itop==0) THEN
            ! No (low) cloud
            emu_flag = 4
            RETURN
        ELSEIF (psrf-pfull(imax)>pmin) THEN
            ! Too high cloud top height
            emu_flag = 5
            RETURN
        ELSEIF (low_iwp>0.25*low_lwp) THEN
            ! No mixed-phase clouds
            emu_flag = 6
            RETURN
        ENDIF

        ! No columns with significant amounts of cloud ice - radiative impacts
        !   - Large ice particles have weak impact on radiative fluxes, for example,
        !     Heymsfield et al. (J. Appl. Meteor., 42, 1369-1390, 2003) found that the visible
        !     optical depth is 2 when IWP=60 g/m2 in mid-latitudes ( tau=0.068*lwp^0.83).
        !     This means transparent ice cloud.
        IF (iwp>60e-3) THEN
            emu_flag = 7
            RETURN
        ENDIF

        ! ******* Emulator inputs for the valid column *******

        ind_lo = MIN(imax+2,klev)
        ind_hi = imax-2
        ! 1) q_inv (g/kg): inversion strength of total water (vapor and cloud) mass mixing ratio
        ! ECHAM: max-min values of the total water for levels from two layers below cloud up to two levels above cloud
        ! IFS: here just water vapor
        q_inv = 1e3*(MAXVAL(qv(ind_hi:ind_lo))-MINVAL(qv(ind_hi:ind_lo)))
        !
        ! 2) tpot_inv (K): inversion strength of (liquid water) potential temperature
        ! ECHAM: max-min values of potential temperature (the same levels as for q_inv)
        tpot_inv = MAXVAL(tpot(ind_hi:ind_lo))-MINVAL(tpot(ind_hi:ind_lo))
        !
        ! 3) tpot_pbl (K): liquid water potential temperature in the boundary layer
        ! ECHAM: min value of potential temperature (the same levels as for q_inv)
        tpot_pbl = MINVAL(tpot(ind_hi:ind_lo))
        !
        ! 4) blh (hPa): inversion/cloud top/planetary boundary layer height described as pressure
        ! ECHAM: pressure difference from surface up to cloud top (top interface of the grid cell)
		! IFS: pressure difference from surface up to cloud to
        pbl = 0.01*(psrf-pfull(imax))
        !
        ! 5) max_lwc (g/kg): maximum liquid water mixing ratio
        ! ECHAM: LWP (g/m2) integrated from the surface up to the cloud top
        ! IFS: IFS: account for cloud coverage (can be zero below cloud) and rain water
        lwc_max = 1e3*(lwc(imax)+rwc(imax))/cc(imax)
        !
        ! 6) CDNC (#/cm3): Cloud droplet number concentration (cdnc) for Level 3 microphysics
        ! ECHAM: averaged over the cloud
        ! IFS: CDNC at the cloud top
        cdnc_max = cdnc(imax)
        !
        ! 7) LWP_rad (g/m2): LWP above the cloud for calcuating radiative fluxes
        lwp_rad=MAX(0.0,(lwp-low_lwp)*1e3)
        !
        ! 8) LHF (W/m2): latent het flux
        lhf = -lhf_ifs ! positive up

        ! ******* Post-filtering *******
        IF (.TRUE.) THEN
            ! No testing
        ELSEIF (q_inv<0.2 .OR. q_inv>10.0) THEN
            ! 0.2 g/kg < q_inv < 10 g/kg
            emu_flag=6
        ELSEIF (tpot_inv<1.0 .OR. tpot_inv>15.0) THEN
            ! 1 K < theta_inv < 15 K
            emu_flag=7
        ELSEIF (tpot_pbl<268.0 .OR. tpot_pbl>305.0) THEN
            ! 268 K < theta < 305 K
            emu_flag=8
        ELSEIF (pbl<25.0) THEN
            ! 25 hPa < BLH
            emu_flag=9
        ELSEIF (lwc_max>1.0) THEN
            ! lwc_max<1 g/kg
            emu_flag=10
        ELSEIF (cdnc_max<10.0) THEN
            ! CDNC > 10 1/mg
            emu_flag=11
        ELSEIF (lwp_rad>300.0) THEN
            ! LWP above <300 g/m2
            emu_flag=12
        ELSEIF (lhf<-50.0 .OR. lhf>400.0) THEN
            ! -50 W/ms < LHF < 400 W/m2
            emu_flag=13
        ENDIF
        !
    END SUBROUTINE test_column_ifs


END MODULE me_emulator
