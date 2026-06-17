!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulator parameters from data files,
! determine emulator input parameters based on column properties,
! and predict updraft velocities.
!
!   Tomi Raatikainen (FMI) 20.3.2026
!
! Current emulator (H25): emulator developments by Noora Hyttinen, 2026.
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
    CHARACTER(LEN=100) :: fname = 'gp_ws_H26.emu'
    ! 2) Emulator input vector
    REAL, SAVE, PUBLIC :: input_vec(8)

    ! Emulator-specific settings
    integer, parameter ::  dp = selected_real_kind(13,300)
    REAL(dp), PARAMETER :: grav = 9.81
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
        ! Local variables for Noora's emulator
        INTEGER :: klevp1
        REAL :: phalf(klev+1)
        !
        ! Emulator variables
        ! a) Dimensions and indices
        klevp1=klev+1
        ! b) Calculate half levels
        phalf(klevp1)=psrf ! Surface
        phalf(2:klev)=0.5*(pfull(2:klev)+pfull(1:klev-1)) ! Between levels
        phalf(1)=0.1*phalf(2) ! Not known, but not relevant
        !
        ! Run emulator
        CALL emulator(klev, klevp1, pfull, phalf, lwc, iwc, qv, cc,shf, cos_mu, tpot, emu_get_inputs)
        !
    END FUNCTION emu_get_inputs


    !>
    !! This subroutine calculates vertical interacls from ground to certain level
    !! used to calculate lwp inside boundary layer
    SUBROUTINE calculate_lwp(x,klev,cloud_top,zdpg,cfrac,integral) 
        INTEGER, INTENT(IN)    :: cloud_top, klev
        REAL(dp),INTENT(IN)    :: x(klev),zdpg(klev),cfrac(klev)
        REAL(dp),INTENT(OUT)   :: integral
        INTEGER                :: jk
        integral = 0
        DO jk = cloud_top,klev
            IF (cfrac(jk) > 0.1) THEN ! Divide by cloud fraction if there is > 10% cloud in that level
                integral  = integral  + x(jk) *zdpg(jk)/cfrac(jk)
            ELSE
                integral  = integral  + x(jk) *zdpg(jk)
            END IF
        END DO
    END SUBROUTINE calculate_lwp

    SUBROUTINE find_cloud2(lev,indhi,x,cloud_top,cloud_base)
        ! Input: lev   = number of layers
        !        indhi = index of the highest layer considered 
        !                (e.g., the 700 hPa level)
        !        x(lev) = cloud liquid water content [g/kg]
        ! Output: cloud_top  = full-level index for the lowermost layer of the lowest cloud 
        !         cloud_base = full-level index for the uppermost layer of the lowest cloud 
        INTEGER, INTENT(IN)   :: lev, indhi
        REAL(dp), INTENT(in)  :: x(lev)
        INTEGER, INTENT(out)  :: cloud_top, cloud_base
        INTEGER :: ilev

        ! Initialize to undefined
        cloud_top=-999
        cloud_base=-999
        DO ilev=lev,indhi,-1
            IF (x(ilev) >= 0.01) THEN
                cloud_base=ilev
                EXIT
            END IF
        ENDDO 
        IF (cloud_base >= indhi) THEN
            DO ilev=cloud_base,indhi,-1
                IF (x(ilev) < 0.01) EXIT
            ENDDO 
            cloud_top=MAX(indhi,ilev+1)
        END IF

    END SUBROUTINE find_cloud2

    !>
    !! This subroutine computes mask where emulator is applied for each timestep.
    !! Returns emulator inputs
    SUBROUTINE emulator( klev, klevp1, &
                      pfull, phalf, pxlm1,pxim1,pqm1,cfrac,shf, cosmu, tpot, pemu_mask)

        ! Input for the emulator 
        INTEGER, INTENT(IN) :: &
          klevp1,         & ! number of half levels (klevp1=klev+1)
          klev              ! number of vertical points, highest number is the surface
        REAL(dp), INTENT(in) :: &
          pfull(klev),    & ! Full-level pressure [Pa]
          phalf(klevp1),  & ! Half-level pressure [Pa]
          pxlm1(klev),    & ! cloud liquid water content [kg/kg]
          pxim1(klev),    & ! cloud ice content [kg/kg]
          pqm1(klev),     & ! specific humidity [kg/kg]
          cfrac(klev),    & ! cloud fraction in levels [0,1]
          shf,            & ! surface sensible heat flux [W/m2]
          cosmu,          & ! solar zenith angle [angle or cos of angle?]
          tpot(klev)        ! potential temperature [K]
        INTEGER, INTENT(OUT) :: pemu_mask    ! emulator mask
        ! Local variables
        INTEGER :: &
          zfull700, &            ! index for nearest level at 700hpa level.
          cloud_top, &           ! index for top of the lowest cloud layer
          cloud_base, &          ! index for base of lowest cloud layer
          inv_ind_hi,inv_ind_lo  ! index's for inversion calculations
        REAL(dp) :: &
          aps, &                 ! surface pressure [Pa]
          zlwp700,ziwp700, &     ! LWP and IWP inside 700hpa layer
          zlwp,ziwp, &           ! Total lwp and iwp
          pdpg(klev), &          ! pressure  difference
          zx(klev),ztotalq(klev) ! zx is just temporary vector,ztotalq is total water

        ! Variables used as the emulator input, the units are ones that are needed for the emulator
        REAL(dp) :: ptpot_inv    ! potential temperature inversion strengt [K]
        REAL(dp) :: ptpot_pbl    ! potential temperature at pbl [K]
        REAL(dp) :: ph2o_inv     ! h2o inversion [kg/kg]
        REAL(dp) :: ppbl_h       ! pbl height [m]

        ! Emulator output
        input_vec(:)=0.0
        pemu_mask = 0

        ! 3. Eliminate all points where there is fog, check if there is cloud water on lowest layer and excluce these points
        IF(pxlm1(klev)*1000 .GE. 0.01) THEN        ! First level water content > 0.01 g/kg
            pemu_mask = 1
            RETURN
        END IF

        ! Emulator is only trained for low level clouds (clouds below 700hpa level). 
        ! 1. step find index for level 700hpa both for half levels and full levels
        zx =pfull(:)
        WHERE(zx .GE. 70000) zx = 999
        zfull700 = MINLOC(MERGE(0,1,zx == 999),DIM=1)

        ! Eliminate points where there is no low level clouds. There is cloud if there is more
        ! than 0.01 g/kg of cloud water below 700hPa
        IF(COUNT(pxlm1(zfull700:klev)*1000 .GE. 0.01) < 1) THEN
            pemu_mask = 2
            RETURN
        END IF

        ! 4. step, locate lowest cloud
        CALL find_cloud2(klev,zfull700,pxlm1(:)*1000,cloud_top,cloud_base)

        IF (.TRUE. .AND. pxlm1(cloud_top-1)*1000 .GE. 0.01) THEN
            ! Bug fix: subroutine find_cloud2 finds the low cloud below 700 hPa, but does not
            ! check if the cloud continues above 700 hPa (the level above)
            pemu_mask = 3
        ENDIF

        ! 5. step calculate LWP inside lowest layer layer, this is used to identify that most of the cloud water in the column is in low level cloud
        ! In-cloud LWP, using cloud fraction of each level (OpenIFS has cloud fraction for each level)
        ! Output zlwp and ziwp are in g/m2
        pdpg(1:klev) = (phalf(2:klevp1) - phalf(1:klev))/grav ! For LWP calculations
        CALL calculate_lwp(pxlm1(:)*1000,klev,cloud_top,pdpg(:),cfrac(:),zlwp700)
        CALL calculate_lwp(pxim1(:)*1000,klev,cloud_top,pdpg(:),cfrac(:),ziwp700)
        IF (.NOT. ( 0.1*zlwp700  > ziwp700)) THEN
            ! Too much ice below 700 hPa
            pemu_mask = 4
            RETURN
        ENDIF
        CALL calculate_lwp(pxlm1(:)*1000,klev,1,pdpg(:),cfrac(:),zlwp)
        CALL calculate_lwp(pxim1(:)*1000,klev,1,pdpg(:),cfrac(:),ziwp)
        IF( .NOT. (zlwp700 > ((zlwp+ziwp)*0.5))) THEN
            ! Clouds above
            pemu_mask = 5
            RETURN
        ENDIF

        ! Acceptable cloud conditions, calculate emulator inputs
        inv_ind_hi = cloud_top-2
        ztotalq = pxlm1(:)+pqm1(:) ! [kg/kg]
        inv_ind_lo = MIN(cloud_base+2,klev)
        ptpot_inv = MAXVAL(tpot(inv_ind_hi:inv_ind_lo))-MINVAL(tpot(inv_ind_hi:inv_ind_lo))
        ph2o_inv = MAXVAL(ztotalq(inv_ind_hi:inv_ind_lo))-MINVAL(ztotalq(inv_ind_hi:inv_ind_lo))
        ptpot_pbl = MINVAL(tpot(inv_ind_hi:inv_ind_lo))
        aps = phalf(klevp1)
        ppbl_h = 44307.69396*((aps/101325)**0.190284-(phalf(cloud_top)/101325)**0.190284) ! in m, positive value

        ! Set emulator input
        input_vec(1)=-shf               ! Sensible heat flux W/m2, positive downwards in OIFS, positive upwards in UCLALES and emulator
        input_vec(2)=cosmu              ! cos of solar zenith angle
        input_vec(3)=aps                ! Surface pressure Pa
        input_vec(4)=ph2o_inv           ! Delta qt kg/kg
        input_vec(5)=ptpot_pbl          ! Theta PBL (cb) K
        input_vec(6)=MAX(ptpot_inv,2.43457) ! Delta theta K, minimum value 2.43K is the lowest value in the training data
        input_vec(7)=zlwp700*1e-3       ! In-cloud LWP kg/m2
        input_vec(8)=ppbl_h             ! Height of PBL (cloud top) m

    END SUBROUTINE emulator

END MODULE me_emulator
