      interface
         function grad(func,v,fixed) result(g)
         use iso_fortran_env, only : r64=>real64
         include "iface/func_iface_inc.f90"
         real(r64) v(:)
         logical, optional :: &
            fixed(:)
         real(r64), allocatable :: &
            g(:)
         end
      end interface
