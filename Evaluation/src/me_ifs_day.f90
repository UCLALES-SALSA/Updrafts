!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read updraft velocity emulators and predict updraft velocities for
! one daily IFS netCDF output
!
!   Tomi Raatikainen (FMI) 12.2.2025
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE me_ifs_day
    USE me_emulator, ONLY : init_emulators, emu_get_inputs, emu_predict, input_vec
    use netcdf
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: IFS_test_day

    ! IFS data from one output file
    CHARACTER(LEN=100) :: fname_in = '../data/training_data_100.nc'
    INTEGER :: nt,nlat,nlon,nlev
    REAL, ALLOCATABLE :: lwc(:,:,:,:),rwc(:,:,:,:),iwc(:,:,:,:),qv(:,:,:,:), &
         press(:,:,:,:),theta(:,:,:,:),cdnc(:,:,:,:),lwp(:,:,:),iwp(:,:,:),cos_mu(:,:,:), &
         lsm(:,:,:),cc(:,:,:,:),pblh(:,:,:),lhfm(:,:,:),shfi(:,:,:), &
         lat(:),lon(:),lev(:),times(:)
    LOGICAL, ALLOCATABLE :: mask(:,:,:)
    ! Default output
    CHARACTER(LEN=100) :: fname_out = 'test_emu_day.dat'

CONTAINS

    SUBROUTINE IFS_test_day(input,output)
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
    END SUBROUTINE IFS_test_day

    ! From sample_data.f90
    SUBROUTINE open_data()
        INTEGER :: i, iret, ncid0, varid
        INTEGER :: timedimid,latdimid,londimid,levdimid
        LOGICAL, SAVE :: first_call=.TRUE.
        !
        ! Open
        iret = nf90_open(trim(fname_in), NF90_NOWRITE, ncid0)
        if (iret /= nf90_noerr) STOP 'Opening failed!'
        !
        ! Dimensions
        iret = nf90_inq_dimid(ncid0,"time",timedimid)
        if (iret /= nf90_noerr) STOP 'Dimension (time) failed!'
        iret = nf90_inq_dimid(ncid0,"lat",latdimid)
        if (iret /= nf90_noerr) STOP 'Dimension (lat) failed!'
        iret = nf90_inq_dimid(ncid0,"lon",londimid)
        if (iret /= nf90_noerr) STOP 'Dimension (lon) failed!'
        iret = nf90_inq_dimid(ncid0,"lev",levdimid)
        if (iret /= nf90_noerr) STOP 'Dimension (lev) failed!'
        !
        IF (first_call) THEN
           ! The first data file - read dimensions and allocate data
           !
           ! Size
           iret = nf90_inquire_dimension(ncid0,timedimid,len=nt)
           iret = nf90_inquire_dimension(ncid0,latdimid,len=nlat)
           iret = nf90_inquire_dimension(ncid0,londimid,len=nlon)
           iret = nf90_inquire_dimension(ncid0,levdimid,len=nlev)
           !
           ! Values
           ALLOCATE(lat(nlat),lon(nlon),lev(nlev),times(nt))
           iret = nf90_inq_varid(ncid0, "lat", varid)
           iret = nf90_get_var(ncid0,varid,lat)
           iret = nf90_inq_varid(ncid0, "lon", varid)
           iret = nf90_get_var(ncid0,varid,lon)
           iret = nf90_inq_varid(ncid0, "lev", varid)
           iret = nf90_get_var(ncid0,varid,lev)
           !
           ! Allocate data: note that dimensions are in the opposite order compared with Python!
           ALLOCATE(lwc(nlon,nlat,nlev,nt),rwc(nlon,nlat,nlev,nt),iwc(nlon,nlat,nlev,nt), &
                qv(nlon,nlat,nlev,nt),press(nlon,nlat,nlev+1,nt),theta(nlon,nlat,nlev,nt), &
                cdnc(nlon,nlat,nlev,nt),lwp(nlon,nlat,nt),iwp(nlon,nlat,nt), &
                cos_mu(nlon,nlat,nt),lsm(nlon,nlat,nt),cc(nlon,nlat,nlev,nt),pblh(nlon,nlat,nt), &
                lhfm(nlon,nlat,nt+1),shfi(nlon,nlat,nt))
           ALLOCATE(mask(nlon,nlat,nt))
           !
           ! Latent heat flux not available for the first time values
           mask(:,:,:) = .TRUE.
           mask(:,:,1) = .FALSE.
           !
           ! Data allocated
           first_call=.FALSE.
        ELSE
           ! Must have the same dimensions
           iret = nf90_inquire_dimension(ncid0,timedimid,len = i)
           if (i /= nt) STOP "Dimension (nt) changes!"
           iret = nf90_inquire_dimension(ncid0,latdimid,len = i)
           if (i /= nlat) STOP "Dimension (nlat) changes!"
           iret = nf90_inquire_dimension(ncid0,londimid,len = i)
           if (i /= nlon) STOP "Dimension (nlon) changes!"
           iret = nf90_inquire_dimension(ncid0,levdimid,len = i)
           if (i /= nlev) STOP "Dimension (nlev) changes!"
           !
           ! All values usable
           mask(:,:,:) = .TRUE.
           !
           ! Latent heat flux is integrated, so to calculate the mean we need the
           ! last value from the previous time step
           lhfm(:,:,1)=lhfm(:,:,nt+1)
        ENDIF
        !
        ! Current time
        iret = nf90_inq_varid(ncid0, "time", varid)
        iret = nf90_get_var(ncid0,varid,times)
        ! 
        ! clwc('time', 'lev', 'lat', 'lon')=lwc: liquid water content (kg/kg)
        iret = nf90_inq_varid(ncid0, "clwc", varid)
        iret = nf90_get_var(ncid0,varid,lwc)
        !
        ! crwc('time', 'lev', 'lat', 'lon')=rwc: rain water content (kg/kg) 
        iret = nf90_inq_varid(ncid0, "crwc", varid)
        iret = nf90_get_var(ncid0,varid,rwc)
        !
        ! ciwc('time', 'lev', 'lat', 'lon')=iwc: ice water content (kg/kg)
        iret = nf90_inq_varid(ncid0, "ciwc", varid)
        iret = nf90_get_var(ncid0,varid,iwc)
        !
        ! q('time', 'lev', 'lat', 'lon')=qv: Specific humidity (kg/kg)
        iret = nf90_inq_varid(ncid0, "q", varid)
        iret = nf90_get_var(ncid0,varid,qv)
        !
        !pt('time', 'lev', 'lat', 'lon')=tpot: Potential temperature (K)
        iret = nf90_inq_varid(ncid0, "pt", varid)
        iret = nf90_get_var(ncid0,varid,theta)
        !
        !pres('time', 'lev', 'lat', 'lon')=pres: Pressure (Pa)
        ! Pressure has one extra level: surface pressure
        iret = nf90_inq_varid(ncid0, "pres", varid)
        iret = nf90_get_var(ncid0,varid,press)
        !
        !param14.216.192('time', 'lev', 'lat', 'lon')=cdnc: CDNC (1/cm3?)
        iret = nf90_inq_varid(ncid0, "param14.216.192", varid)
        iret = nf90_get_var(ncid0,varid,cdnc)
        !
        !var78('time', 'lat', 'lon')=lwp: Total column liquid water (kg/m2)
        iret = nf90_inq_varid(ncid0, "var78", varid)
        iret = nf90_get_var(ncid0,varid,lwp)
        !
        !var79('time', 'lat', 'lon')=iwp: Total column ice water (kg/m2)
        iret = nf90_inq_varid(ncid0, "var79", varid)
        iret = nf90_get_var(ncid0,varid,iwp)
        !
        !param26.216.192('time', 'lev', 'lat', 'lon')=cos_mu: Cosine of solar zenith angle (-)
        ! Note: level-dependent cos_mu => use the data from the lowest level
        iret = nf90_inq_varid(ncid0, "param26.216.192", varid)
        iret = nf90_get_var(ncid0,varid,cos_mu,start=(/1,1,nlev,1/),count=(/nlon,nlat,1,nt/))
        !
        !var172('time', 'lat', 'lon')=lsm: Land-sea mask (-)
        iret = nf90_inq_varid(ncid0, "var172", varid)
        iret = nf90_get_var(ncid0,varid,lsm)
        !
        !cc('time', 'lev', 'lat', 'lon')=cc: Cloud cover (-)
        iret = nf90_inq_varid(ncid0, "cc", varid)
        iret = nf90_get_var(ncid0,varid,cc)
        !
        !var159('time', 'lat', 'lon')=pblh: Boundary layer height (m)
        iret = nf90_inq_varid(ncid0, "var159", varid)
        iret = nf90_get_var(ncid0,varid,pblh)
        !
        !var147('time', 'lat', 'lon')=lhfm: Surface latent heat flux (J/m2)
        ! Note: lhfm(:,:,1) contains the last value from the previous data file (day)
        iret = nf90_inq_varid(ncid0, "var147", varid)
        iret = nf90_get_var(ncid0,varid,lhfm(:,:,2:))
        !
        !var231('time', 'lat', 'lon')=shfi: Instantaneous surface sensible heat flux (W/m2)
        iret = nf90_inq_varid(ncid0, "var231", varid)
        iret = nf90_get_var(ncid0,varid,shfi)

        ! Ignored:
        !param13.212.192('time', 'lev', 'lat', 'lon')=nait: number dry aitken (1/g)
        !param16.212.192('time', 'lev', 'lat', 'lon')=nacc: number dry accumulation (1/g)
        !param22.212.192('time', 'lev', 'lat', 'lon')=ncoa: number dry coarse (1/g)

        iret = nf90_close(ncid0)

    END SUBROUTINE open_data


    SUBROUTINE main_loop()
        INTEGER :: i, j, k, n, iout
        INTEGER :: n_good, n_bad, n_land, emu_flag, flags(20)
        REAL :: pred
        LOGICAL :: fail
        !
        n_good = 0
        n_bad = 0
        n_land = 0
        flags(:) = 0
        !
        ! Output dump
        OPEN(newunit=iout,FILE=fname_out)
        !
        DO i=1,nt
         DO j=1,nlat
          do k=1,nlon
            IF (mask(k,j,i) .AND. lsm(k,j,i)<0.1 .AND. lwp(k,j,i)>0.0 .AND. theta(k,j,nlev,i)>265.0) THEN
                ! Determine inputs (input_vec) for the current column
                emu_flag = emu_get_inputs(nlev, press(k,j,2:,i), theta(k,j,:,i), lwc(k,j,:,i), &
                        rwc(k,j,:,i), iwc(k,j,:,i), cc(k,j,:,i), qv(k,j,:,i), cdnc(k,j,:,i), lwp(k,j,i), &
                        iwp(k,j,i), press(k,j,1,i), lhfm(k,j,i), shfi(k,j,i), lsm(k,j,i), cos_mu(k,j,i))
                !
                IF (emu_flag==0) THEN
                    ! Valid column, so try to predict
                    CALL emu_predict(input_vec,pred,fail)
                    IF (fail) THEN
                        ! Inputs not with the range of training data
                        n_bad = n_bad+1
                        flags(15) = flags(15) + 1
                    ELSE
                        n_good = n_good+1
                    ENDIF
                else
                    n_bad = n_bad+1
                ENDIF
            ELSE
                n_land = n_land+1
            ENDIF
          ENDDO
         ENDDO
        END DO
        CLOSE(13)
        CLOSE(14)
        !
        ! Report
        n=n_land+n_good+n_bad
        WRITE(*,'(A7,I7,F7.2)') "Land: ",n_land,n_land/(0.01*N)
        WRITE(*,'(A7,I7,F7.2)') "Good: ",n_good,n_good/(0.01*N)
        WRITE(*,'(A7,I7,F7.2)') "Bad: ", n_bad,n_bad/(0.01*N)
        DO i=1,SIZE(flags)
            IF (flags(i)>0) WRITE(*,'(I7,I7,F7.2)') i,flags(i),flags(i)/(0.01*N)
        ENDDO
        !
    END SUBROUTINE main_loop

END MODULE me_ifs_day