! ****************************************
! Original source: GPF from https://github.com/ots22/gpf
! Modified 11.2.2025 by Tomi Raatikainen (FMI)
! ****************************************
module me_optim
  use me_gp, only: DenseGP, nlog_lik
  implicit none
  include 'nlopt.f'

contains


  subroutine log_lik_optim(nhypers, gp, lbounds, ubounds, max_niter, ftol_rel, maxf, niter)
    TYPE(DenseGP), intent(in) :: gp
    integer, intent(in) :: nhypers, max_niter
    real, intent(in) :: lbounds(nhypers), ubounds(nhypers), ftol_rel
    real, intent(out) :: maxf
    integer, intent(out) :: niter 
    real :: hypers(nhypers)
    integer(kind=8) :: opt
    integer nnu, ires

    ! Hyperparameters: just theta or both nu and theta
    if (nhypers==size(gp%theta)) then
        hypers = gp%theta
    else
        nnu = size(gp%nu)
        hypers(1:nnu) = gp%nu
        hypers(nnu+1:) = gp%theta
    endif

    call nlo_create(opt, NLOPT_LN_BOBYQA, nhypers)
    call nlo_set_lower_bounds(ires, opt, lbounds)
    call     check_error_code(ires)
    call nlo_set_upper_bounds(ires, opt, ubounds)
    call     check_error_code(ires)
    call nlo_set_max_objective(ires, opt, nlog_lik, gp)
    call     check_error_code(ires)
    call nlo_set_ftol_rel(ires, opt, ftol_rel)
    call     check_error_code(ires)
    call nlo_set_maxeval(ires, opt, max_niter)
    call     check_error_code(ires)
    call nlo_optimize(ires, opt, hypers, maxf)
    call     check_error_code(ires)
    call nlo_get_numevals(niter, opt)
    call nlo_destroy(opt)

  end subroutine log_lik_optim

  subroutine check_error_code(ires)
    integer ires
    if (ires<0) then
       print *, "NLopt failed with error code ", ires
       stop 1
    end if   
  end subroutine check_error_code

end module me_optim
