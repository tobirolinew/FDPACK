      interface
         subroutine diag_hess(func,v,eigvals,eigvects,fixed)
         use iso_fortran_env, only : r64=>real64
         include "iface/func_iface_inc.f90"
         real(r64) :: v(:)
         real(r64), allocatable :: &
            h(:,:),eigvals(:),eigvects(:,:)
         logical, optional :: &
            fixed(:)
         end
      end interface
