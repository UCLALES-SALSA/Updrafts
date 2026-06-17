! ****************************************
! Original source: GPF from https://github.com/ots22/gpf
! GP (m_gp.f90 and m_gp_dense.f90)
!
! Modified 11.2.2025 by Tomi Raatikainen (FMI)
!
! A Gaussian process of full rank.  See [1, section 2], which gives an
! introduction to Gaussian processes.
! 
! [1] J. Quinonero and C. Rasmussen. Analysis of Some Methods for
! Reduced Rank Gaussian Process Regression, in Switching and Learning
! in Feedback Systems: European Summer School on Multi-Agent Control,
! Maynooth, Ireland, September 8-10, 2003, Revised Lectures and
! Selected Papers, Springer, 2005

module me_gp
  use me_util, only : pi, solve
  implicit none
  private
  PUBLIC DenseGP, update_matrices, nlog_lik, nnu_required, ntheta_required, write_out, predict

  type :: DenseGP
     ! The noise hyperparameter(s). This is passed to the noise model
     ! (`noise_model'), which determines its precise meaning.
     real, dimension(:), allocatable :: nu 
     ! covariance hyperparameters
     ! meaning depends on covariance function `covariance'
     real, dimension(:), allocatable :: theta
     ! inputs
     real, dimension(:,:), allocatable :: x
     ! observations
     real, dimension(:), allocatable :: t
     ! precomputed product, used in the prediction
     real, dimension(:), allocatable :: invCt
     ! LU decomposition of the covariance matrix, used in optimization
     real, dimension(:,:), allocatable, PRIVATE :: LU
     !
     ! Predictions require some information about the training data:
     ! minimum and maximum values (the range of training data),  and
     ! the means and standard deviations for normalization
     REAL, ALLOCATABLE :: min_x(:), max_x(:), mean_x(:), std_x(:)
     REAL :: mean_t = 0.0, std_t = 1.0
  end type DenseGP

  interface DenseGP
     module procedure make_DenseGP
     module procedure read_DenseGP
  end interface DenseGP

  contains

    ! ****************************************
    ! GP (m_gp.f90)
  
    ! A routine with the required interface for the NLopt library, for
    ! maximizing the log-likelihood.
    subroutine nlog_lik(val, n, hypers, grad, need_gradient, gp)
        TYPE(DenseGP) :: gp
        integer :: n, need_gradient
        real :: val, hypers(n)
        real, intent(inout) :: grad(n)
        integer :: nnu, ntheta
        real :: logprior

        if (need_gradient/=0) then
            print *, "optimization requires gradients of hyperparameters"
            stop 1
        end if

        nnu = size(gp%nu)
        ntheta = size(gp%theta)

        if (ntheta==n) THEN
            ! hypers=[theta]
            gp%theta = hypers
        ELSEIF (nnu+ntheta==n) then
            ! hypers=[nu,theta]
            gp%nu = hypers(1:nnu)
            gp%theta = hypers(nnu+1:nnu+ntheta)
        ELSE
            print *, "size of optimization variables does not match the number of hyperparameters"
            stop 1
        end if

        ! Update 
        call update_matrices(gp)

        ! The old log-likelihood:
        !logprior = - log(2. * sqrt(gp%theta(1)))
        !do i = 2,ntheta
        !    logprior = logprior - log(1. + gp%theta(i)**2)
        !end do

        ! The new log-likelihood
        logprior = -0.5*size(gp%t)*log(2.0*pi)

        val = log_lik(gp) + logprior

        call output_params(gp,val,logprior)
    end subroutine nlog_lik

    ! Helper routine for nlog_lik: print out the current hyperparameters
    ! (associated with `gp') and their corresponding log-likelihood (`val').
    subroutine output_params(gp,val,logprior)
        TYPE(DenseGP) :: gp
        real :: val,logprior
        integer, save :: u, niter=0
        character(len=100), SAVE :: conv

        if (niter==0) THEN
            open(newunit=u, file="LOG_LIK_OPTIM")
            WRITE(conv, FMT="(A11,I2,A10,I2,A17)") "(A9,I5,A10,",SIZE(gp%nu), &
                "F9.6,A10,",SIZE(gp%theta),"F12.7,A26,2F12.5)"
        ENDIF

        niter=niter+1

        write (u,conv) "Iteration",niter,": noise = ", gp%nu, ", theta = ", & 
             gp%theta, ", log(prior,likelihood) = ", logprior,val

    end subroutine output_params



    ! ****************************************
    ! Dense GP (m_gp_dense.f90)
    subroutine alloc_DenseGP(gp, n, ntheta, dims)
        type(DenseGP), intent(inout) :: gp
        ! n: number of observations
        ! ntheta: number of covariance hyperparameters
        ! dims: dimension of the inputs
        integer, intent(in) :: n, ntheta, dims

        allocate(gp%nu(nnu_required(dims)))
        allocate(gp%theta(ntheta))
        allocate(gp%x(n, dims))
        allocate(gp%t(n))
        allocate(gp%invCt(n))
        allocate(gp%LU(n,n))

        ALLOCATE(gp%min_x(dims), gp%max_x(dims))
        ALLOCATE(gp%mean_x(dims), gp%std_x(dims))

    end subroutine alloc_DenseGP

    function make_DenseGP(nu, theta, x, t, normalize) result(gp)
        type(DenseGP) :: gp
        ! noise hyperparameters
        real, dimension(:), intent(in) :: nu
        ! covariance hyperparameters
        real, dimension(:), intent(in) :: theta
        ! training input coordinates
        real, dimension(:,:), intent(in) :: x
        ! training outputs
        real, dimension(:), intent(in) :: t
        ! normalize training data
        logical, intent(in) :: normalize

        integer :: n, d
        n = size(t)
        d = size(x,2)

        if (size(theta) /= ntheta_required(d)) then
            print *, "size of theta does not match number of hyperparameters required by the covariance function"
            stop 1
        end if

        if (size(nu) /= nnu_required(d)) then
            print *, "size of nu (noise params) does not match number required by the noise model"
            stop 1
        end if

        call alloc_DenseGP(gp, n, size(theta), d)

        gp%nu = nu
        gp%theta = theta
        gp%t = t
        gp%x = x

        if (normalize) then
            call NormalizeTraining(gp)
        else
            gp%min_x=-HUGE(1.0); gp%max_x=HUGE(1.0)
            gp%mean_x=0.0; gp%std_x=1.0
            gp%mean_t=0.0; gp%std_t=1.0
        endif

        call update_matrices(gp)

    end function make_DenseGP

    subroutine write_out(this, filename)
        TYPE(DenseGP), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: u ! unit number for output

        open(newunit=u, file=filename)

        write (u,'(A)') "DenseGPE"
        write (u,'(I10)') size(this%t), size(this%theta), size(this%nu), size(this%x,2)
        write (u,'(es24.15)') this%nu, this%theta, this%x, this%invCt
        write (u,'(es24.15)') this%min_x, this%max_x
        write (u,'(es24.15)') this%mean_x, this%std_x, this%mean_t, this%std_t
     
        close(u)
    end subroutine write_out

    function read_DenseGP(filename) result(gp)
        character(len=*), intent(in) :: filename
        type(DenseGP) :: gp
        integer n, ntheta, nnu, d, u
        character(len=12) :: label

        open(newunit=u, file=filename)
        read (u,'(A)') label
        if (trim(label) /= "DenseGPE") then
           print *, "Incompatible data file"
           stop 1
        end if
        read (u,'(I10)') n, ntheta, nnu, d
        call alloc_DenseGP(gp, n, ntheta, d)
        read (u,'(es24.15)') gp%nu, gp%theta, gp%x, gp%invCt
        read (u,'(es24.15)') gp%min_x, gp%max_x
        read (u,'(es24.15)') gp%mean_x, gp%std_x, gp%mean_t, gp%std_t
        close(u)

        ! This is a new variable for log-likelihood
        gp%LU = 0.

    end function read_DenseGP


    ! Helper routine to update the internal state.  Called when an
    ! observation or hyperparameter changes and the covariance matrix
    ! must be recomputed.
    subroutine update_matrices(this)
        TYPE(DenseGP), intent(inout) :: this
        integer :: i,j,n
        ! always a small amount of noise to stabilize the inversion
        real, parameter :: noise_stab = 1e-9
        real :: C(size(this%t),size(this%t))

        n = size(this%t)

        ! Q in ref [1] (just below eq. (7) therein)
        do i=1,n
           do j=1,i-1
              C(i,j) = cov(this%x(i,:), this%x(j,:), this%theta)
              C(j,i) = C(i,j)
           end do
           C(i,i) = cov(this%x(i,:), this%x(i,:), this%theta) + noise(this%nu) + noise_stab
        end do

        this%invCt = solve(C, this%t, this%LU)
    end subroutine update_matrices

     ! Make a prediction of the underlying function value at
     ! coordinate `xnew'.
     function predict(this, xnew)
        real :: predict
        TYPE(DenseGP), intent(in) :: this
        real, dimension(:), intent(in) :: xnew
        integer :: i
        real :: k(size(this%t)), x(size(xnew))

        ! Normalize
         x = (xnew - this%mean_x)/this%std_x

        do i=1,size(this%t)
           k(i) = cov(x,this%x(i,:),this%theta)
        end do
        ! Predictive mean from eq. (6) in ref [1]. Here, k corresponds to
        ! k^{*} in the reference, and t to y
        predict = dot_product(k, this%invCt)*this%std_t+this%mean_t ! Unnormalize
     end function predict

    ! The log likelihood of the hyperparameters
    function log_lik(this)
        real :: log_lik
        TYPE(DenseGP), intent(in) :: this
        ! First line of eq. (8) in ref. [1]
        log_lik = -0.5 * (logdet(this%LU) + dot_product(this%t, this%invCt))
     end function log_lik

    ! Calculates log determinant based on the LU decomposition (originally from util.f90)
    function logdet(LU)
        real :: logdet
        real, dimension(:,:), intent(in) :: LU
        integer :: i
        logdet = 0.0
        do i=1,size(LU,1)
            logdet = logdet + log(abs(LU(i,i)))
        end do
    end function logdet


    ! ****************************************
    ! Noise model (noise_value_only.f90)

    pure function nnu_required(dims)
        integer nnu_required
        integer, intent(in) :: dims
        ! single noise level to be applied to the target value (and zero to the derivatives)
        nnu_required = 1
    end function nnu_required

    pure function noise(params)
        real :: noise
        real, intent(in) :: params(:)
        noise = params(1)
    end function noise


    ! ****************************************
    ! Covariance function (cov_linsqexp.f90)

    pure function ntheta_required(dims)
        integer :: ntheta_required
        integer, intent(in) :: dims
        ! the scale parameter plus one `r' parameter for each dimension
        ntheta_required = dims+1
    end function ntheta_required

    pure function cov(x,y,hypers)
        real :: cov
        real, intent(in), dimension(:) :: x, y, hypers
        real :: scale, r(size(x,1))
        scale = hypers(1)
        r(:) = hypers(2:)
        !cov = scale * exp(-0.5*sum((x-y)**2/r**2)) + 10 + 10 * sum(x * y) ! Nordling et al., (2024)
        !cov = scale + exp(-0.5*sum((x-y)**2/r**2)) ! Ahola et al. (2022)
        cov = scale * exp(-0.5*sum((x-y)**2/r**2)) ! New
    end function cov


    ! ****************************************
    ! Normalize training data and return normalization coefficients (mean and
    ! standard deviation) for further use. Also return the range (min and max)
    ! of the training data inputs.

    SUBROUTINE NormalizeTraining(this)
        TYPE(DenseGP), intent(inout) :: this
        ! Local
        INTEGER :: i, n, m
        !
        ! Calculate minimum and maximum values of the unnormalized inputs
        this%min_x=MINVAL(this%x,DIM=1)
        this%max_x=MAXVAL(this%x,DIM=1)
        !
        ! Calculate means and standard deviations, and normalize
        n=SIZE(this%x,1)
        m=SIZE(this%x,2)
        DO i=1,m
            this%mean_x(i) = mean(this%x(:,i),n)
            this%std_x(i) = std(this%x(:,i),this%mean_x(i),n)
            this%x(:,i) = (this%x(:,i) - this%mean_x(i))/this%std_x(i)
        ENDDO
        this%mean_t = mean(this%t,n)
        this%std_t = std(this%t,this%mean_t,n)
        this%t(:) = (this%t(:) - this%mean_t)/this%std_t
        !
        CONTAINS
            FUNCTION mean(x,dmn) RESULT(res)
                INTEGER :: dmn
                REAL :: x(dmn)
                REAL :: res
                res = SUM(x)/dmn
            END FUNCTION mean
            !
            FUNCTION std(x,meanx,dmn) RESULT(res)
                INTEGER :: dmn
                REAL :: x(dmn)
                REAL :: meanx
                REAL :: res
                res = SQRT(SUM((x - meanx)**2)/dmn)
            END FUNCTION std
    END SUBROUTINE NormalizeTraining

end module me_gp
