!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Emulator training and testing.
! Modified 11.2.2025 by Tomi Raatikainen (FMI)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE me_training
    USE me_gp, ONLY : DenseGP, update_matrices, nnu_required, ntheta_required, predict, write_out
    use me_optim, ONLY : log_lik_optim
    IMPLICIT NONE
    !
    PRIVATE
    PUBLIC leave_one_out, leave_some_out, test_emulator
    !
    ! Emulator training data to be read from a text file
    REAL, ALLOCATABLE :: ex2(:,:), t2(:)
    ! Emulator
    TYPE(DenseGP), ALLOCATABLE :: gp
CONTAINS


    ! ************ Tests ************

    SUBROUTINE read_training(fname,n,m)
        ! Function for reading training data text files into variables ex2(:,:) and t2(:)
        CHARACTER(*), INTENT(IN) :: fname ! Data file name
        INTEGER, INTENT(OUT) :: n, m ! Dimensions
        ! Emulator inputs and outputs
        INTEGER :: uu, i
        CHARACTER :: com
        !
        ! Input file
        open(newunit=uu,file=fname,ACTION="READ")
        !
        ! New files have the dimensions (rows, columns) in the first line
        read (uu,*) n,m
        ! The second line may contain column headers after a comment character
        READ (uu,'(A1)') com
        IF (com/='!' .AND. com/='#') BACKSPACE(uu)
        !
        ! Number of variables
        m=m-1
        ! Allocate emulator training data
        IF (ALLOCATED(ex2)) DEALLOCATE(ex2,t2)
        ALLOCATE(ex2(n,m),t2(n))
        ! Read all
        DO i=1,n
            read (uu,*) ex2(i,1:m), t2(i)
        ENDDO
        close(uu)
        WRITE(*,*) 'Data loaded, dimensions',n,m
        if (n .le. 10) STOP 'Not enough data!'
    END SUBROUTINE read_training


    SUBROUTINE test_emulator(dataname,emuname,outname,outemu,iseed)
        ! Compare existing or new emulator against training data
        ! Inputs
        CHARACTER(*), INTENT(IN) :: dataname,emuname ! Training data and emulator files
        CHARACTER(*), INTENT(IN) :: outname ! Output file name
        CHARACTER(*), INTENT(IN) :: outemu  ! Output emulator name
        INTEGER, INTENT(IN), OPTIONAL :: iseed ! Optional random number generator seed
        ! Local variables
        REAL, ALLOCATABLE :: emuInputVec(:)
        INTEGER :: i, n, m, iout
        REAL :: pred, rmse, aad
        !
        ! Read training data: ex2(:,:) and t2(:)
        CALL read_training(dataname,n,m)
        !
        ! Allocate data
        ALLOCATE(emuInputVec(m))
        !
        IF (LEN_TRIM(emuname)>1) THEN
            ! Read the emulator
            IF (allocated(gp)) DEALLOCATE(gp)
            allocate(gp, source = DenseGP(emuname))
        ELSE
            ! Train emulator
            !CALL train_emulator_rand(0, n, m, ex2, t2) ! Use the default as an initial guess
            CALL train_emulator_rand(30, n, m, ex2, t2, iseed) ! Random initial values
        ENDIF
        !
        ! Output
        iout=6 ! Negative if file opened, otherwise standard output
        IF (LEN_TRIM(outname)>1) OPEN(newunit=iout,FILE=outname,ACTION='WRITE')
        !
        ! Predict
        WRITE(*,*) ' '
        WRITE(*,*) 'Test with the training data set'
        rmse=0.
        aad=0.
        DO i=1,n
            emuInputVec(:)=ex2(i,1:m)
            ! Predict
            pred = predict(gp,emuInputVec)
            ! Print
            WRITE(iout,*) i, t2(i), pred, t2(i)-pred
            ! RMSE and AAD
            rmse=rmse+(t2(i)-pred)**2
            aad=aad+abs(t2(i)-pred)
        ENDDO
        IF (iout<0) CLOSE(iout)
        WRITE(*,'(A7,F10.5)') 'RMSE: ',SQRT(rmse/FLOAT(n))
        WRITE(*,'(A7,F10.5)') 'AAD:  ',aad/FLOAT(n)
        !
        ! Save emulator only if it is trained here (emuname='')
        IF (LEN_TRIM(outemu)>1 .AND. LEN_TRIM(emuname)==0) THEN
            WRITE(*,*) ' '
            WRITE(*,*)'Saving parameters to '//TRIM(outemu)
            WRITE(*,'(A7,F10.5)') '  nu',  gp%nu
            WRITE(*,'(A7,15F10.5,/)') 'theta', gp%theta
            call write_out(gp,outemu)
        ENDIF
        !
        ! Clean
        DEALLOCATE(emuInputVec)
    END SUBROUTINE test_emulator


    SUBROUTINE leave_some_out(frac,dataname,emuname,outname,iseed)
        ! Leave-some-out tests where
        ! Inputs
        REAL, INTENT(IN) :: frac ! Fraction of the test data set from the total number of points
        CHARACTER(*), INTENT(IN) :: dataname ! Training data file
        CHARACTER(*), INTENT(IN) :: emuname ! Input emulator name (optional)
        CHARACTER(*), INTENT(IN) :: outname ! Output file name
        INTEGER, INTENT(IN), OPTIONAL :: iseed ! Optional random number generator seed
        ! Local variables
        REAL, ALLOCATABLE :: ex1(:,:), t1(:), emuInputVec(:), th0(:), nu0(:)
        REAL :: log_lik_val
        INTEGER :: i, j, k, n_test, nn, n, m, iout, niter
        INTEGER, ALLOCATABLE :: ind(:)
        REAL, ALLOCATABLE :: pred(:)
        !
        ! Read training data: ex2(:,:) and t2(:)
        CALL read_training(dataname,n,m)
        !
        IF (LEN_TRIM(emuname)>1) THEN
            ! Read previously trained emulator - will be used as the initial value
            IF (allocated(gp)) DEALLOCATE(gp)
            allocate(gp, source = DenseGP(emuname))
        ELSE
            ! Find the best fit for all the data - will be used as the initial value
            ! Training with random initial values
            CALL train_emulator_rand(10, n, m, ex2, t2, iseed)
        ENDIF
        !
        ! Save the initial hyperparameter values
        ALLOCATE(th0(size(gp%theta)),nu0(size(gp%nu)))
        th0=gp%theta
        nu0=gp%nu
        !
        ! Allocate index and outputs
        ALLOCATE(ind(n),pred(n))
        !
        ! Index vector with values ranging from 1 to x so that each index is included n_test times
        n_test = ceiling(frac*n) ! The last test data set can be smaller than n_test
        CALL generate_rand_ind(n_test,n,ind)
        !
        ! Output
        iout=6 ! Negative if file opened, otherwise standard output
        IF (LEN_TRIM(outname)>1) OPEN(newunit=iout,FILE=outname,ACTION='WRITE')
        !
        WRITE(*,*)' '
        WRITE(*,*)'Training emulators for leave-n-out test...'
        WRITE(*,*)'ind  pnts niter     fval     nu       th0      th1      th2      ...'
        !
        ! Leave-n_test-out tests with n-n_test observations
        ALLOCATE(ex1(n-n_test,m),t1(n-n_test),emuInputVec(m))
        DO i=1,MAXVAL(ind)
            ! Number of training data values
            nn=COUNT(ind/=i) ! All other points
            ! Training data
            IF (nn/=n-n_test) THEN
                DEALLOCATE(ex1,t1)
                ALLOCATE(ex1(nn,m),t1(nn))
            ENDIF
            k=1
            DO j=1,n
                IF (ind(j)/=i) THEN
                     ex1(k,1:m)=ex2(j,1:m)
                     t1(k)=t2(j)
                     k=k+1
                ENDIF
            ENDDO
            !
            ! Reset the emulator
            CALL allocate_default(nn, m, ex1, t1)
            ! Initial values
            gp%theta=th0
            gp%nu=nu0
            !
            ! Emulator training
            CALL train_emulator(niter, log_lik_val)
            !
            WRITE(UNIT=*,FMT='(I4,2I6,F12.2,15F9.4)') i, nn, niter, log_lik_val, gp%nu, gp%theta
            !
            ! Predictions
            DO j=1,n
                IF (ind(j)==i) THEN
                    ! Emulator input (row j)
                    emuInputVec(:)=ex2(j,1:m)
                    ! Predict
                    pred(j) = predict(gp,emuInputVec)
                ENDIF
            ENDDO
        ENDDO

        ! Write output (the LES output, emulator prediction and the difference
        WRITE(*,*)' '
        WRITE(*,*)'Predictions'
        IF (iout>0) WRITE(*,*)'  i       obs      pred       err ind' ! Screen prints
        DO i=1,n
            WRITE(iout,'(I4,3F10.4,I4)') i, t2(i), pred(i), t2(i)-pred(i), ind(i)
        ENDDo
        ! RMSE, AAD
        WRITE(*,*) 'RMSE: ',SQRT(SUM((t2(:)-pred(:))**2)/FLOAT(n))
        WRITE(*,*) 'AAD:  ',SUM(ABS((t2(:)-pred(:)))/FLOAT(n))

        ! All done
        IF (iout<0) CLOSE(iout)
        DEALLOCATE(ex1, t1, emuInputVec, th0, nu0, ind, pred, gp)
    END SUBROUTINE leave_some_out
    !
    ! Generates index vector ivec(n) so that it contains random indexes from 1 to ceil(n/m).
    ! Each index is included m times, except that the last is inluded mod(n,m) times.
    SUBROUTINE generate_rand_ind(m,n,ivec)
        INTEGER, INTENT(IN) :: m, n
        INTEGER, INTENT(OUT) :: ivec(n)
        REAL :: rands(n)
        INTEGER :: i, k(1), ind, nn
        !
        ! Random numbers from 0 to 1
        call random_number( rands )
        !
        ! Indexes (1,2,..)
        ind=1 ! Index (1,2,..)
        nn=0 ! Number of index values
        DO i=1,n
            ! Find the lowest value (0-1)
            k = MINLOC(rands)
            ! Current index
            ivec(k(1))=ind
            rands(k(1))=2.0 ! Exclude
            nn=nn+1
            ! No more than m
            IF (m==nn) THEN
                ind=ind+1
                nn=0
            ENDIF
        ENDDO
    END SUBROUTINE generate_rand_ind

    SUBROUTINE leave_one_out(dataname,emuname,outname,iseed)
        ! Leave-one-out tests
        ! Inputs
        CHARACTER(*), INTENT(IN) :: dataname ! Training data file
        CHARACTER(*), INTENT(IN) :: emuname ! Input emulator name (optional)
        CHARACTER(*), INTENT(IN) :: outname ! Output file name
        INTEGER, INTENT(IN), OPTIONAL :: iseed ! Optional random number generator seed
        ! Local variables
        REAL, ALLOCATABLE :: ex1(:,:), t1(:), emuInputVec(:), th0(:), nu0(:)
        REAL :: pred, log_lik_val
        INTEGER :: i, n, m, iout, niter
        !
        ! Read training data: ex2(:,:) and t2(:)
        CALL read_training(dataname,n,m)
        !
        IF (LEN_TRIM(emuname)>1) THEN
            ! Read previously trained emulator - will be used as the initial value
            IF (allocated(gp)) DEALLOCATE(gp)
            allocate(gp, source = DenseGP(emuname))
        ELSE
            ! Find the best fit for all the data - will be used as the initial value
            ! Training with random initial values
            CALL train_emulator_rand(10, n, m, ex2, t2, iseed)
        ENDIF
        !
        ! Save the best values
        ALLOCATE(th0(size(gp%theta)),nu0(size(gp%nu)))
        th0=gp%theta
        nu0=gp%nu
        !
        !
        ! Output
        iout=0 ! Negative if file opened, otherwise standard output
        IF (LEN_TRIM(outname)>1) OPEN(newunit=iout,FILE=outname,ACTION='WRITE')
        !
        ! Leave-one-out tests with n-1 observations
        ALLOCATE(ex1(n-1,m),t1(n-1),emuInputVec(m)) 
        DO i=1,n
            ! Training data (without row i)
            IF (i>1) THEN
                ex1(1:i-1,1:m)=ex2(1:i-1,1:m)
                t1(1:i-1)=t2(1:i-1)
            ENDIF
            IF (i<n) THEN
                ex1(i:n-1,1:m)=ex2(i+1:n,1:m)
                t1(i:n-1)=t2(i+1:n)
            ENDIF
            ! Emulator input (row i)
            emuInputVec(:)=ex2(i,1:m)
            !
            ! Reset the emulator
            CALL allocate_default(n-1, m, ex1, t1)
            ! Initial values
            gp%theta=th0
            gp%nu=nu0
            !
            ! Emulator training
            CALL train_emulator(niter, log_lik_val)
            !
            ! Predict
            pred = predict(gp,emuInputVec)
            !
            ! Write output (the LES output, emulator prediction and the difference, and 
            ! the total number if iterations and the final function value; hyperparameters if ouput file)
            WRITE(UNIT=*,FMT=*) i, t2(i), pred, t2(i)-pred, niter, log_lik_val
            IF (iout<0) WRITE(iout,*) i, t2(i), pred, t2(i)-pred, niter, log_lik_val, gp%nu, gp%theta
        ENDDO
        ! All done
        IF (iout<0) CLOSE(iout)
        DEALLOCATE(ex1, t1, emuInputVec, th0, nu0, gp)
    END SUBROUTINE leave_one_out


    SUBROUTINE train_emulator_rand(nrepeats, n, m, x, t, iseed)
        ! Emulator training with m different initial values drawn from random distributions
        ! Inputs
        integer, INTENT(IN) :: nrepeats ! Number of  repeats
        integer, INTENT(IN) :: n, m ! Dimensions of the inputs
        real, INTENT(IN) :: x(n,m), t(n) ! Training data
        INTEGER, INTENT(IN), OPTIONAL :: iseed ! Optional random number generator seed
        ! Emulator parameters
        INTEGER :: nnu, ntheta, i, j, k(1), niter
        REAL, ALLOCATABLE :: nu(:,:), theta(:,:)
        REAL :: fval(nrepeats), diff, log_lik_val
        integer, allocatable, dimension(:) :: seed
        !
        ! Generate repeatable pseudorandom numbers
        IF (PRESENT(iseed)) THEN
            call random_seed(size=j)
            allocate (seed(j))
            seed = iseed * (/ (i, i = 1, j) /) + 5
            call random_seed(put=seed)
        ENDIF
        !
        ! Allocate emulator - the same x and t for all repeats
        CALL allocate_default(n, m, x, t)
        !
        IF (nrepeats<=0) THEN
            ! Single fit based on the default parameters
            CALL train_emulator(niter, log_lik_val)
            WRITE(*,*)' '
            WRITE(*,*)'Emulator training with default initial values...'
            WRITE(*,*)'Niter       fval      nu     th0     th1     th2 ...'
            WRITE(*,'(I6,F11.3,F8.4,15F8.3)') niter, log_lik_val, gp%nu(:), gp%theta(:)
            RETURN
        ENDIF
        !
        ! Allocate nu and theta for all repeats
        nnu = SIZE(gp%nu) ! Noise parameters
        ntheta = SIZE(gp%theta) ! Other parameters
        ALLOCATE(nu(nrepeats,nnu),theta(nrepeats,ntheta))
        j=1 ! Current best
        WRITE(*,*)' '
        WRITE(*,*)'Emulator training with random initial values...'
        WRITE(*,*)'  Niter       fval      nu     th0     th1     th2 ...'
        DO i=1,nrepeats
            ! Generate random theta
            call random_number( theta(i,:) ) ! 0-1
            !theta(i,:) = 2.0*theta(i,:)+0.1 ! 0.1-2.1
            !theta(i,:) = 4.0*theta(i,:)+0.1 ! 0.1-4.1
            theta(i,:) = 10.0**( 2.*theta(i,:)-1.0 ) ! log10(theta) is evenly distributed between -1 and 1
            !theta(i,:) = get_random_tri(0.1,10.0,1.0,ntheta) ! triangular, 0.1-10
            !theta(i,:) = get_random_tri(0.1,5.0,1.0,ntheta) ! triangular, 0.1-5
            ! Generate random nu
            call random_number( nu(i,:) ) ! 0-1
            !
            ! Set the initial values
            gp%nu=nu(i,:)
            gp%theta=theta(i,:)
            !
            ! Train
            CALL train_emulator(niter, log_lik_val)
            WRITE(*,'(I2,I6,F11.3,F8.4,15F8.3)') i, niter, log_lik_val, gp%nu(:), gp%theta(:)
            !
            ! Is this the best?
            IF (log_lik_val>fval(j)) j=i
            !
            ! Save parameters
            fval(i) = log_lik_val
            theta(i,:) = gp%theta(:)
            nu(i,:) = gp%nu(:)
        ENDDO
        !
        ! Set current emulator to the best
        gp%nu=nu(j,:)
        gp%theta=theta(j,:)
        call update_matrices(gp)
        !
        ! Return if less than two repeats
        IF (nrepeats<2) RETURN
        !
        ! Print sorted results
        WRITE(*,*)' '
        WRITE(*,*)'The best results'
        diff = MAXVAL(fval)-MINVAL(fval)
        DO j=1,MIN(10,nrepeats)
            ! Current maximum
            k = MAXLOC(fval)
            i =  k(1)
            ! Print
            WRITE(UNIT=*,FMT='(2I4,F11.3,F8.4,15F8.3)') j, i, fval(i), nu(i,:), theta(i,:)
            ! Change to new minimum
            fval(i) = fval(i) - 2*diff
        ENDDO
        !
    END SUBROUTINE train_emulator_rand
    !
    ! Get random number vector from a triangular distribution
    !   Range from a to b and peak at c
    FUNCTION get_random_tri(a,b,c,n) RESULT(rnd)
        REAL :: a,b,c !  a<c<b
        INTEGER :: n
        REAL :: rnd(n) 
        INTEGER :: i
        !
        call random_number( rnd ) ! Evenly distibuted values from 0 to 1
        DO i=1,n
            IF ( rnd(i)<(c-a)/(b-a) ) THEN
                rnd(i) = a + SQRT(rnd(i)*(b-a)*(c-a))
            ELSE
                rnd(i) = b - SQRT((1.-rnd(i))*(b-a)*(b-c))
            END IF
        END DO
        !
    END FUNCTIOn get_random_tri


    SUBROUTINE train_emulator(niter, log_lik_val)
        ! Emulator training
        INTEGER, INTENT(OUT) :: niter
        REAL, INTENT(OUT) :: log_lik_val
        ! Local
        INTEGER :: nnu, ntheta
        ! Optimization parameter
        INTEGER :: max_iter = 10000
        real :: ftol = 1.0e-7
        REAL, SAVE, ALLOCATABLE :: lbounds(:),ubounds(:)

        ! Emulator must exists
        IF (.not.allocated(gp)) STOP 'Emulator not allocated!'

        ! Bounds for hyperparameter optimization: allocate if not yet allocated
        IF (allocated(lbounds)) THEN
            ! Use the current values
        ELSEIF (.FALSE.) THEN
            ! Noise parameter is fixed
            nnu = SIZE(gp%nu)
            ntheta = SIZE(gp%theta)
            ALLOCATE(lbounds(ntheta),ubounds(ntheta))
            lbounds(:) = 0.01
            ubounds(:) = 100.0
        ELSE
            ! Noise parameter is the first hyperparameter
            nnu = SIZE(gp%nu)
            ntheta = SIZE(gp%theta)
            ALLOCATE(lbounds(nnu+ntheta),ubounds(nnu+ntheta))
            lbounds(:) = 0.01
            ubounds(:) = 100.0
            ! Noise from 1e-4 to 1
            lbounds(1) = 1e-4
            ubounds(1) = 1.0
        ENDIF

        ! Optimization
        call log_lik_optim(SIZE(lbounds), gp, lbounds, ubounds, max_iter, ftol, log_lik_val, niter)

    END SUBROUTINE train_emulator


    SUBROUTINE allocate_default(n, m, x, t)
        ! Inputs
        integer, INTENT(IN) :: n, m ! Number of training points and dimension of the input
        real, INTENT(IN) :: x(n,m), t(n) ! Training data
        ! Emulator parameters
        INTEGER :: nnu, ntheta
        REAL, ALLOCATABLE :: nu(:),theta(:)
        LOGICAL :: normalize
        !
        !  Defaults
        nnu = nnu_required(m) ! Noise parameters
        ntheta = ntheta_required(m) ! Other parameters
        ALLOCATE(nu(nnu),theta(ntheta))
        nu(:) = 0.001
        theta(:) = 1.
        !
        ! Normalize training data
        normalize = .TRUE.
        !
        ! Reset emulator
        IF (allocated(gp)) DEALLOCATE(gp)
        allocate(gp, source=DenseGP(nu, theta, x, t, normalize))
        !
        ! Clean
        DEALLOCATE(nu,theta)
    END SUBROUTINE allocate_default

END MODULE me_training
