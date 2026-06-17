!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulators and predict updraft velocities for
! the IFS sample dataset either based on the netCDF dataset or the
! previously calculated emulator inputs
!
!   Tomi Raatikainen (FMI) 12.2.2025
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE me_ifs_sample
    USE me_emulator, ONLY : init_emulators, emu_get_inputs, emu_predict, input_vec
    use netcdf
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: IFS_sample_full

    ! IFS sample dataset
    !CHARACTER(LEN=100) :: fname_in = '../data/sample_dataset.nc' ! Marine columns
    CHARACTER(LEN=100) :: fname_in = '../data/sample_dataset_all.nc' ! Land and marine columns
    INTEGER :: nt,nlev
    REAL,ALLOCATABLE :: cdnc(:,:),theta(:,:),press(:,:),qv(:,:),lwc(:,:),iwc(:,:), &
                rwc(:,:),cc(:,:),lwp(:),iwp(:),lsm(:),pblh_ifs(:),lhfm(:),shfi(:), &
                cos_mu(:),p_surf(:),lev(:),lat(:),lon(:),times(:)
    ! Default output
    CHARACTER(LEN=100) :: fname_out = 'test_emu_sample.dat'

CONTAINS

    ! Compute emulator inputs and predictions
    SUBROUTINE IFS_sample_full(input,output)
        character (len=*), intent(in), optional :: input, output
        ! Input and output files
        if (present(input)) fname_in = input ! NetCDF input
        if (present(output)) fname_out = output ! Text output file
        ! Read the input data
        CALL open_data()
        ! Setup the default emulator
        CALL init_emulators()
        ! Predictions
        CALL main_loop()
    END SUBROUTINE IFS_sample_full

    ! From process_sample.f90
    SUBROUTINE open_data()
        INTEGER :: iret, ncid, varid
        !
        ! Open
        iret = nf90_open(trim(fname_in), NF90_NOWRITE, ncid)
        if (iret /= nf90_noerr) STOP 'Opening failed!'
        !
        ! Dimensions
        iret = nf90_inq_dimid(ncid,"time",varid)
        if (iret /= nf90_noerr) STOP 'Dimensions (time) failed!'
        iret = nf90_inquire_dimension(ncid,varid,len = nt)
        iret = nf90_inq_dimid(ncid,"lev",varid)
        if (iret /= nf90_noerr) STOP 'Dimensions (lev) failed!'
        iret = nf90_inquire_dimension(ncid,varid,len = nlev)
        !
        ! Allocate data: note that dimensions are in the opposite order compared with Python!
        ALLOCATE(cdnc(nlev,nt),theta(nlev,nt),press(nlev,nt),qv(nlev,nt),lwc(nlev,nt),iwc(nlev,nt), &
                rwc(nlev,nt),cc(nlev,nt),lwp(nt),iwp(nt),lsm(nt),pblh_ifs(nt),lhfm(nt),shfi(nt),&
                cos_mu(nt),p_surf(nt),lev(nlev),lat(nt),lon(nt),times(nt))
        !
        ! Time, coordinates and level indices
        iret = nf90_inq_varid(ncid,"time",varid)
        iret = nf90_get_var(ncid,varid,times)
        iret = nf90_inq_varid(ncid,"lev",varid)
        iret = nf90_get_var(ncid,varid,lev)
        iret = nf90_inq_varid(ncid,"lat",varid)
        iret = nf90_get_var(ncid,varid,lat)
        iret = nf90_inq_varid(ncid,"lon",varid)
        iret = nf90_get_var(ncid,varid,lon)
        !
        ! 1D outputs
        ! ========
        ! var78('time', 'lat', 'lon')=lwp: Total column liquid water (kg/m2)
        iret = nf90_inq_varid(ncid,"lwp",varid)
        iret = nf90_get_var(ncid,varid,lwp)
        ! var79('time', 'lat', 'lon')=iwp: Total column ice water (kg/m2)
        iret = nf90_inq_varid(ncid,"iwp",varid)
        iret = nf90_get_var(ncid,varid,iwp)
        ! var172('time', 'lat', 'lon')=lsm: Land-sea mask (-)
        iret = nf90_inq_varid(ncid,"lsm",varid)
        iret = nf90_get_var(ncid,varid,lsm)
        ! var159('time', 'lat', 'lon')=pblh: Boundary layer height (m)
        iret = nf90_inq_varid(ncid,"pblh",varid)
        iret = nf90_get_var(ncid,varid,pblh_ifs)
        ! var147('time', 'lat', 'lon')=lhfm: Mean surface latent heat flux (W/m2)
        iret = nf90_inq_varid(ncid,"lhfm",varid)
        iret = nf90_get_var(ncid,varid,lhfm)
        ! var231('time', 'lat', 'lon')=shfi: Instantaneous surface sensible heat flux (W/m2)
        iret = nf90_inq_varid(ncid,"shfi",varid)
        iret = nf90_get_var(ncid,varid,shfi)
        ! param26.216.192('time', 'lev', 'lat', 'lon')=cos_mu: Cosine of solar zenith angle (-)
        iret = nf90_inq_varid(ncid,"cos_mu",varid)
        iret = nf90_get_var(ncid,varid,cos_mu)
        ! Surface pressure taken from the first level of pressure array
        iret = nf90_inq_varid(ncid,"p_surf",varid)
        iret = nf90_get_var(ncid,varid,p_surf)
        !
        ! 2D outputs
        ! ========
        ! param14.216.192('time', 'lev', 'lat', 'lon')=cdnc: CDNC (1/cm3)
        iret = nf90_inq_varid(ncid,"cdnc",varid)
        iret = nf90_get_var(ncid,varid,cdnc)
        ! pt('time', 'lev', 'lat_2', 'lon_2')=tpot: Potential temperature (K)
        iret = nf90_inq_varid(ncid,"theta",varid)
        iret = nf90_get_var(ncid,varid,theta)
        ! pres('time', 'lev', 'lat_2', 'lon_2')=pres: Pressure (Pa)
        iret = nf90_inq_varid(ncid,"press",varid)
        iret = nf90_get_var(ncid,varid,press)
        ! q('time', 'lev', 'lat', 'lon')=qv: Specific humidity (kg/kg)
        iret = nf90_inq_varid(ncid,"qv",varid)
        iret = nf90_get_var(ncid,varid,qv)
        ! clwc('time', 'lev', 'lat', 'lon')=lwc: liquid water content (kg/kg)
        iret = nf90_inq_varid(ncid,"lwc",varid)
        iret = nf90_get_var(ncid,varid,lwc)
        ! ciwc('time', 'lev', 'lat', 'lon')=iwc: ice water content (kg/kg)
        iret = nf90_inq_varid(ncid,"iwc",varid)
        iret = nf90_get_var(ncid,varid,iwc)
        ! crwc('time', 'lev', 'lat', 'lon')=rwc: rain water content (kg/kg)
        iret = nf90_inq_varid(ncid,"rwc",varid)
        iret = nf90_get_var(ncid,varid,rwc)
        ! cc('time', 'lev', 'lat', 'lon')=cc: Cloud cover (-)
        iret = nf90_inq_varid(ncid,"cc",varid)
        iret = nf90_get_var(ncid,varid,cc)
        !
        iret = nf90_close(ncid)
    END SUBROUTINE open_data

    SUBROUTINE main_loop()
        INTEGER :: i, iout, emu_flag
        INTEGER :: n, n_good, n_bad, flags(20)=0
        REAL :: pred
        LOGICAL :: fail
        !
        n_good = 0
        n_bad = 0
        flags(:) = 0
        !
        ! Output emulator predictions
        OPEN(newunit=iout,FILE=fname_out)
        !
        DO i=1,nt
            ! Determine inputs (input_vec) for the current column
            emu_flag = emu_get_inputs(nlev, press(:,i), theta(:,i), lwc(:,i), rwc(:,i), iwc(:,i), cc(:,i), &
                                 qv(:,i), cdnc(:,i), lwp(i), iwp(i), p_surf(i), lhfm(i), shfi(i), lsm(i), cos_mu(i))
            !
            IF (emu_flag==0) THEN
                ! Valid column, so try to predict
                CALL emu_predict(input_vec,pred,fail)
                !
                IF (fail) THEN
                    ! Inputs not within the range of training data
                    n_bad = n_bad+1
                    emu_flag = 15
                    flags(emu_flag) = flags(emu_flag)+1
                ELSE
                    n_good = n_good+1
                ENDIF
            ELSE
                ! Failed
                n_bad = n_bad+1
                flags(emu_flag)=flags(emu_flag)+1
                pred=-999.
            ENDIF
            WRITE(iout,*) pred, emu_flag
        END DO
        !
        ! Report
        n=n_good+n_bad
        WRITE(*,'(A7,I7,F7.2)') "Good: ",n_good,n_good/(0.01*N)
        WRITE(*,'(A7,I7,F7.2)') "Bad: ", n_bad,n_bad/(0.01*N)
        DO i=1,SIZE(flags)
            IF (flags(i)>0) WRITE(*,'(I7,I7,F7.2)') i,flags(i),flags(i)/(0.01*N)
        ENDDO
        !
    END SUBROUTINE main_loop

END MODULE me_ifs_sample