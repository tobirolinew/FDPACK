!***********************************************************************
      module mod_fdbase  ! base module for finite differences
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      implicit none
      type fdconfig_
         character(7), public :: &
            !finite-difference scheme
            !('central' or 'forward'):
               fdscheme='central'
         real(r64), public, allocatable :: &
            !array of finite-difference offsets:
            ! (i,i,1): gradient step for variable i
            ! (i,i,2): diagonal Hessian step for variable i
            ! (i,j,1): i-direction step for entry (i,j)
            ! (i,j,2): j-direction step for entry (i,j)
               fdoffsets(:,:,:)
      end type
      type(fdconfig_) :: &
         fdconfig
      private
      public :: fdconfig, get_num_grad, get_num_hess, fwgrad,&
                cdgrad, fwhess, cdhess, diag_num_hess
      !$omp threadprivate(fdconfig)
      contains
      !*****************************************************************
         function get_num_grad(func,v,fixed) result(g)
      !*****************************************************************
      !
      !  This function computes the numerical gradient vector of a
      !  multivariate function by finite forward/central differences.
      !  It uses the global vector "offsets".
      !  An optional mask vector allows fixing selected variables.
      !
      !  Returns:
      !    g(:) - gradient of "func" at "v"
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
            !type for local variables:
         type lview_
            logical, allocatable :: fixed(:)
         end type
         type(lview_) :: &
            !lview_ instance:
               lview
         logical, optional :: &
            !flags for variable fixation:
               fixed(:)
         real(r64) :: &
            !variable vector:
               v(:),&
            !function value at "v"
               f
         real(r64), allocatable :: &
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !gradient of "func" at "v":
               g(:)
         integer :: &
            !iterator:
               i
         !***
         associate(nv=>size(v))
            !initialize the variables in "lview"
            if(present(fixed)) then
               lview%fixed=fixed
            else
               lview%fixed=spread(.false.,1,nv)
            endif

            !check for size incompatibility
            call catch_error(&
                     err=nv.ne.size(lview%fixed).or.&
                         nv.ne.size(fdconfig%fdoffsets,1).or.&
                         nv.ne.size(fdconfig%fdoffsets,2).or.&
                         size(fdconfig%fdoffsets,3).ne.2,&
                     msg='size incompatibility detected.',&
                     proc='get_num_grad')

            !set the current offsets, zeroing fixed variables
            fdoffsets=fdconfig%fdoffsets
            do i=1,nv
               if(.not.lview%fixed(i)) cycle
               fdoffsets(i,i,1)=0.d0
            end do

            !compute the gradient at "v"
            select case(fdconfig%fdscheme)
               case('central')
                  g=cdgrad(func,v,fdoffsets)
               case('forward')
                  f=func(v)
                  g=fwgrad(func,v,fdoffsets,f)
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown fdscheme.',&
                           proc='get_num_grad')
            end select
         end associate
         !***
         end

      !*****************************************************************
         function get_num_hess(func,v,fixed) result(h)
      !*****************************************************************
      !
      !  This function computes the numerical Hessian matrix of a
      !  multivariate function by finite forward/central differences.
      !  It uses the global "offsets" vector.
      !  An optional mask vector allows fixing selected variables.
      !
      !  Returns:
      !    h(:,:) - Hessian of "func" at "v"
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
         logical, optional :: &
            !flags for variable fixation:
               fixed(:)
         real(r64) :: &
            !variable vector:
               v(:)
            !type for local variables:
         type lview_
            real(r64), allocatable :: v(:)
            logical, allocatable :: fixed(:)
         end type
         type(lview_) :: &
            !lview_ instance:
               lview
         real(r64) :: &
            !function value at "v":
               f
         real(r64), allocatable :: &
            !gradient vector at "v":
               g(:),&
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !Hessian matrix at "v":
               h(:,:)
         integer :: &
            !iterator:
               i
         !***
         associate(nv=>size(v))
            !initialize the variables in "lview"
            if(present(fixed)) then
               lview%fixed=fixed
            else
               lview%fixed=spread(.false.,1,nv)
            endif

            !check for size incompatibility
            call catch_error(&
                     err=nv.ne.size(lview%fixed).or.&
                         nv.ne.size(fdconfig%fdoffsets,1).or.&
                         nv.ne.size(fdconfig%fdoffsets,2).or.&
                         size(fdconfig%fdoffsets,3).ne.2,&
                     msg='size incompatibility detected.',&
                     proc='get_num_hess')

            !set the current offsets, zeroing fixed variables
            fdoffsets=fdconfig%fdoffsets
            do i=1,nv
               if(.not.lview%fixed(i)) cycle
               fdoffsets(i,:,:)=0.d0
               fdoffsets(:,i,:)=0.d0
            end do

            !compute the function value at "v"
            f=func(v)

            !compute the Hessian at "v"
            select case(fdconfig%fdscheme)
               case('central')
                  h=cdhess(func,v,fdoffsets,f)
               case('forward')
                  g=fwgrad(func,v,fdoffsets,f)
                  h=fwhess(func,v,fdoffsets,f,g)
               case default
                  call catch_error(&
                           err=.true.,&
                           msg='unknown fdscheme.',&
                           proc='get_num_hess')
            end select
         end associate
         !***
         end

      !*****************************************************************
         subroutine diag_num_hess(func,v,eigvals,eigvects,fixed)
      !*****************************************************************
      !
      !  This subroutine computes the spectral decomposition of the
      !  Hessian matrix for a multivariate function.
      !  The Hessian is determined numerically by calling
      !  "get_num_hess", using the global "offsets" vector.
      !  An optional mask vector allows fixing selected variables.
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
            !Lapack interfaces:
         include 'src/iface/lapack_iface_inc.f90'
         logical, optional :: &
            !flags for variable fixation:
            fixed(:)
         real(r64) :: &
            !variable vector:
               v(:)
         integer :: &
            !error flag for MKL DSYEV:
               info
         real(r64), allocatable :: &
            !workspace for MKL DSYEV:
               work(:),&
            !array of eigenvalues:
               eigvals(:),&
            !array of eigenvectors (and also the Hessian):
               eigvects(:,:)
         !***
         associate(nv => size(v))
            !calculate the Hessian and store it in "eigvects"
            eigvects=get_num_hess(func,v,fixed)

            !initialize "eigvals" and "work"
            eigvals=spread(0.d0,1,nv)
            work=spread(0.d0,1,3*nv-1)

# ifdef USE_LAPACK
            !call the symmetric diagonalizer
            call dsyev('V','U',nv,eigvects,nv,eigvals,work,&
                       size(work),info)
# else
            !report error: USE_LAPACK is not set
            call catch_error(&
                  err=.true., &
                  msg='diag_num_hess requires USE_LAPACK.', &
                  proc='diag_num_hess')
# endif
            !check for errors
            call catch_error(&
                     err=info.ne.0,&
                     msg='DSYEV failed.',&
                     proc='diag_num_hess')
         end associate
         !***
         end

      !*****************************************************************
         function cdgrad(func,v,fdoffsets,f,hdiag) result(g)
      !*****************************************************************
      !
      !  This function computes the gradient of a multivariate function
      !  using finite central differences.
      !
      !  Returns:
      !    g(:)     - gradient of "func" at "v"
      !    hdiag(:) - diagonal of the Hessian at "v" (optional)
      !
      !*****************************************************************
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
         real(r64) :: &
            !variable vector:
               v(:),&
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !function values at perturbed points:
               fp,fm
         real(r64), allocatable :: &
            !temporary variable vector:
               vtmp(:),&
            !gradient of "func" at "v":
               g(:)
         real(r64), optional, allocatable :: &
            !Hessian diagonal:
               hdiag(:)
         real(r64), optional :: &
            !function value at "v"
               f
         integer :: &
            !iterator:
               i
         !***
         associate(nv=>size(v))
            !fill up the temporary variable vector
            vtmp=v

            !initialize "g"
            g=spread(0.d0,1,nv)

            !if present, initialize "hdiag"
            if(present(f).and.present(hdiag)) &
               hdiag=spread(0.d0,1,nv)

            !iterate over the variables
            do i=1,nv
               associate(vi    => v(i),&
                         vti   => vtmp(i),&
                         offi  => fdoffsets(i,i,1),&
                         hoffi => 0.5d0*fdoffsets(i,i,1))
                  !skip entry if the offset is zero
                  if(offi.eq.0.d0) cycle

                  !evaluate the forward step
                  vti=vi+hoffi
                  fp=func(vtmp)

                  !evaluate the backward step
                  vti=vi-hoffi
                  fm=func(vtmp)

                  !save the gradient component
                  g(i)=(fp-fm)/offi

                  !save the Hessian diagonal
                  if(present(f).and.present(hdiag)) &
                     hdiag(i)=(fp-2.d0*f+fm)/hoffi**2

                  !reset the ith variable
                  vti=vi
               end associate
            end do
         end associate
         !***
         end

      !*****************************************************************
         function fwgrad(func,v,fdoffsets,f) result(g)
      !*****************************************************************
      !
      !  This function computes the gradient of a multivariate function
      !  using finite forward differences.
      !
      !  Returns:
      !    g(:)     - gradient of "func" at "v"
      !
      !*****************************************************************
         implicit none
            !objective function
         include 'src/iface/func_iface_inc.f90'
         real(r64) :: &
            !variable vector:
               v(:),&
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !function value at "v"
               f,&
            !function values at the perturbed points:
               fp
         real(r64), allocatable :: &
            !temporary variable vector:
               vtmp(:),&
            !gradient of "func" at "v":
               g(:)
         integer :: &
            !iterator:
               i
         !***
         associate(nv=>size(v))
            !fill up the temporary variable vector
            vtmp=v

            !initialize the gradient vector
            g=spread(0.d0,1,nv)

            !iterate over the variables
            do i=1,nv
               associate(vi   => v(i),&
                         vti  => vtmp(i),&
                         offi => fdoffsets(i,i,1))
                  !skip entry if the offset is zero
                  if(offi.eq.0.d0) cycle

                  !evaluate the forward step
                  vti=vi+offi
                  fp=func(vtmp)

                  !save the gradient component
                  g(i)=(fp-f)/offi

                  !reset the ith variable
                  vti=vi
               end associate
            end do
         end associate
         !***
         end

      !*****************************************************************
         function cdhess(func,v,fdoffsets,f) result(h)
      !*****************************************************************
      !
      !  This function computes the Hessian of a multivariate function
      !  using finite central differences (direct function calls).
      !
      !  Returns:
      !    h(:,:) - Hessian of "func" at "v"
      !
      !*****************************************************************
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
         real(r64) :: &
            !variable vector:
               v(:),&
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !function value at "v":
               f,&
            !perturbed function values:
               fp,fm,fpp,fpm,fmp,fmm
         real(r64), allocatable :: &
            !temporary variable vector:
               vtmp(:),&
            !Hessian of "func" at "v":
               h(:,:)
         integer :: &
            !iterators:
            i, j
         !***
         associate(nv=>size(v))
            !fill up the temporary vector
            vtmp=v

            !initialize the Hessian
            allocate(h(nv,nv),source=0.d0)

            !iterate over the variables
            do i=1,nv
               associate(vi  => v(i),&
                         vti => vtmp(i))
                  !compute the diagonal entry
                  associate(hoffi=>0.5d0*fdoffsets(i,i,2))
                     if(hoffi.ne.0.d0) then
                        vtmp(i)=vi+hoffi
                        fp=func(vtmp)
                        vtmp(i)=vi-hoffi
                        fm=func(vtmp)
                        vtmp(i)=vi
                        h(i,i)=(fp-2.d0*f+fm)/hoffi**2
                     endif
                  end associate

                  !compute the off-diagonal entries
                  do j=1,i-1
                     associate(vj    => v(j),&
                               vtj   => vtmp(j),&
                               hoffi => 0.5d0*fdoffsets(i,j,1),&
                               hoffj => 0.5d0*fdoffsets(i,j,2))
                        !skip the actual entry if fixed
                        if(hoffi.eq.0.d0.or.hoffj.eq.0.d0) cycle

                        !perturbation (+,+)
                        vti=vi+hoffi
                        vtj=vj+hoffj
                        fpp=func(vtmp)

                        !perturbation (+,-)
                        vtj=vj-hoffj
                        fpm=func(vtmp)

                        !perturbation (-,-)
                        vti=vi-hoffi
                        fmm=func(vtmp)

                        !perturbation (-,+)
                        vtj=vj+hoffj
                        fmp=func(vtmp)

                        !save the Hessian entry
                        h(i,j)=0.25d0*(fpp-fpm-fmp+fmm)/hoffi/hoffj
                        h(j,i)=h(i,j)

                        !reset the actual variables
                        vti=vi
                        vtj=vj
                     end associate
                  end do
               end associate
            end do
         end associate
         !***
         end

      !*****************************************************************
         function fwhess(func,v,fdoffsets,f,g) result(h)
      !*****************************************************************
      !
      !  This function computes the Hessian of a multivariate function
      !  using finite forward differences (direct function calls).
      !
      !  Returns:
      !    h(:,:) - Hessian of "func" at "v"
      !
      !*****************************************************************
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
         real(r64) :: &
            !variable vector:
               v(:),&
            !gradient at "v":
               g(:),&
            !finite-difference offsets:
               fdoffsets(:,:,:),&
            !function value at "v":
               f,&
            !perturbed function values:
               fpp
         real(r64), allocatable :: &
            !temporary variable vector:
               vtmp(:),&
            !Hessian of "func" at "v":
               h(:,:)
         integer :: &
            !iterators:
               i,j
         !***
         associate(nv=>size(v))
            !fill up the temporary vector
            vtmp=v

            !initialize the Hessian
            allocate(h(nv,nv),source=0.d0)

            !iterate over the variables
            do i=1,nv
               associate(vi  => v(i),&
                         vti => vtmp(i),&
                         gi  => g(i))
                  !iterate over the lower triangular entries
                  do j=1,i
                     associate(vj   => v(j),&
                               vtj  => vtmp(j),&
                               offi => fdoffsets(i,j,1),&
                               offj => fdoffsets(i,j,2),&
                               gj   => g(j))
                        !skip the actual entry if fixed
                        if(offi.eq.0.d0.or.offj.eq.0.d0) cycle

                        !diagonal entry with identical offsets
                        vti=vti+offi
                        vtj=vtj+offj
                        fpp=func(vtmp)
                        vti=vi
                        vtj=vj
                        h(i,j)=(fpp-f)/offi/offj-gi/offj-gj/offi
                        h(j,i)=h(i,j)
                     end associate
                  end do
               end associate
            end do
         end associate
         end
      end
