! ****************************************
! Original source: GPF from https://github.com/ots22/gpf
! Modified 11.2.2025 by Tomi Raatikainen (FMI)
! ****************************************
module me_util
  implicit none
  real, parameter :: pi = 3.14159265358979323846264338327

contains
  ! Solve a system of linear equations A*x = b for a matrix A and
  ! vectors x and b.  Implemented as a wrapper around LAPACK dgesv.
  function solve(M,b,LU) result(x)
    real, intent(in) :: M(:,:), b(:)
    real, intent(OUT) :: LU(:,:)
    real, dimension(size(b,1)) :: x
    integer :: ipiv(size(b,1)), N, info

    ! check M is square and b conforms
    if (size(M,1) /= size(M,2) .OR. size(M,1) /= size(b,1)) then
       write (*,*) "solve: Array M passed to solve should be square.  Got ", &
            size(M,1), "x", size(M,2)
       error stop
    end if
    N = size(M, 1)

    ! LU and x are overwritten on output of dgesv, so copy
    LU = M
    x = b
    call dgesv(N, 1, LU, N, ipiv, x, N, info)

    ! check for success
    if (info /= 0) then
       write (*,*) "solve: dgesv returned an error code (", info, ")"
       error stop
    end if

  end function solve

end module me_util
