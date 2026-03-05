!***********************************************************************
      program main
!***********************************************************************
!     DESCRIPTION:
!***********************************************************************
!
!     This program demonstrates the full workflow of FDPACK on a
!     two-variable test function:
!
!        f(x,y)=sin(a*x)*exp(-b*y^2)
!
!     Steps covered:
!       [1] noise estimation       (get_noise_table)
!       [2] offset selection       (get_fdoffsets, central + forward)
!       [3] numerical gradient     (get_num_grad)
!       [4] numerical Hessian      (get_num_hess)
!       [5] spectral decomposition (diag_num_hess)
!
!     The 'fixed' mask is set to .false. for all variables.
!     Setting any entry to .true. would freeze that variable
!     (zero gradient/Hessian row/col), useful for constrained
!     problems.
!
!***********************************************************************
      use iso_fortran_env, only : &
         r64=>real64
      use mod_fdbase, only : &
         fdconfig, get_num_grad, get_num_hess,&
         diag_num_hess
      use mod_fdtune, only : &
         get_noise_table, get_fdoffsets
      implicit none
      integer :: &
         ncalls=0
      real(r64), parameter :: &
         a=3.d0,&
         b=2.d0
      real(r64) :: &
         v(2),&
         f0,&
         g_ana(2),&
         h_ana(2,2),&
         g_err(2),&
         h_err(2,2)
      real(r64), allocatable :: &
         noise_table(:,:),&
         fnoises(:),&
         minoffsets(:),&
         maxoffsets(:),&
# ifdef USE_LAPACK         
         eigvals(:),&
         eigvects(:,:),&
# endif
         g(:),&
         h(:,:)
      logical :: &
         fixed(2)
      character(len=16) :: &
         lbl
      integer :: &
         i,j
      !***
      v(1)=0.5d0
      v(2)=0.5d0
      fixed=.false.
      maxoffsets=[1.d0,1.d0]
      !--- analytic reference ---
      f0=func(v)
      g_ana(1)=a*dcos(a*v(1))*dexp(-b*v(2)**2)
      g_ana(2)=-2.d0*b*v(2)*dsin(a*v(1))*dexp(-b*v(2)**2)
      h_ana(1,1)=-a**2*dsin(a*v(1))*dexp(-b*v(2)**2)
      h_ana(2,1)=-2.d0*b*v(2)*a*dcos(a*v(1))*&
                 dexp(-b*v(2)**2)
      h_ana(1,2)=h_ana(2,1)
      h_ana(2,2)=dsin(a*v(1))*dexp(-b*v(2)**2)*&
                 (-2.d0*b+4.d0*b**2*v(2)**2)
      !================================================================
      print'(a)',repeat('=',60)
      print'(a)',' FDPACK example'
      print'(a)',' f(x,y)=sin(a*x) * exp(-b*y^2)'
      print'(a,f4.1,a,f4.1)',&
            '  a=',a,'   b=',b
      print'(a,f6.3,a,f6.3)',&
            '  evaluation point: x=',v(1),'   y=',v(2)
      print'(a,es12.4)',' f(x,y)=',f0
      print'(a)',repeat('=',60)
      !================================================================
      print'(/,a)',' [1] noise estimation'
      print'(a)',repeat('-',60)
      print'(a)',' calling get_noise_table ...'
      noise_table=get_noise_table(func,v,maxoffsets)
      minoffsets=noise_table(:,1)
      fnoises=noise_table(:,2)
      print'(a,2es12.4)',' estimated noise  : ',&
            fnoises(1),fnoises(2)
      print'(a,2es12.4)',' min offsets      : ',&
            minoffsets(1),minoffsets(2)
      print'(a,2es12.4)',' max offsets      : ',&
            maxoffsets(1),maxoffsets(2)
      !================================================================
      print'(/,a)',&
            ' [2a] offset selection  (central, approx=.false.)'
      print'(a)',repeat('-',60)
      fdconfig%fdscheme='central'
      print'(a)',' calling get_fdoffsets ...'
      ncalls=0
      fdconfig%fdoffsets=&
         get_fdoffsets(&
            func,v,fnoises=fnoises,&
            minoffsets=minoffsets,&
            maxoffsets=maxoffsets,&
                       approx=.false.)
      print'(a,i0)',' function calls   : ',ncalls
      print'(a)',' optimal offsets:'
      do i=1,2
         associate(oi=>fdconfig%fdoffsets(i,i,1))
            write(lbl,'(a,i0,a)') 'g(',i,')'
            print'(4x,a16,a,es12.4)',lbl,' ->',oi
         end associate
      end do
      do i=1,2
         do j=1,i
            if(i.eq.j) then
               associate(od=>fdconfig%fdoffsets(i,i,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row/col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',od
               end associate
            else
               associate(&
                     or=>fdconfig%fdoffsets(i,j,1),&
                     oc=>fdconfig%fdoffsets(i,j,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row'
                  print'(4x,a16,a,es12.4)',lbl,' ->',or
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',oc
               end associate
            end if
         end do
      end do
      !================================================================
      print'(/,a)',&
            ' [2b] offset selection  (central, approx=.true.)'
      print'(a)',repeat('-',60)
      fdconfig%fdscheme='central'
      print'(a)',' calling get_fdoffsets ...'
      ncalls=0
      fdconfig%fdoffsets=&
         get_fdoffsets(&
            func,v,fnoises=fnoises,&
            minoffsets=minoffsets,&
            maxoffsets=maxoffsets,&
            approx=.true.)
      print'(a,i0)',' function calls   : ',ncalls
      print'(a)',' optimal offsets:'
      do i=1,2
         associate(oi=>fdconfig%fdoffsets(i,i,1))
            write(lbl,'(a,i0,a)') 'g(',i,')'
            print'(4x,a16,a,es12.4)',lbl,' ->',oi
         end associate
      end do
      do i=1,2
         do j=1,i
            if(i.eq.j) then
               associate(od=>fdconfig%fdoffsets(i,i,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row/col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',od
               end associate
            else
               associate(&
                     or=>fdconfig%fdoffsets(i,j,1),&
                     oc=>fdconfig%fdoffsets(i,j,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row'
                  print'(4x,a16,a,es12.4)',lbl,' ->',or
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',oc
               end associate
            end if
         end do
      end do
      !================================================================
      print'(/,a)',&
            ' [2c] offset selection  (forward, approx=.false.)'
      print'(a)',repeat('-',60)
      fdconfig%fdscheme='forward'
      print'(a)',' calling get_fdoffsets ...'
      ncalls=0
      fdconfig%fdoffsets=&
         get_fdoffsets(&
         func,v,fnoises=fnoises,&
         minoffsets=minoffsets,&
         maxoffsets=maxoffsets,&
         approx=.false.)
      print'(a,i0)',' function calls   : ',ncalls
      print'(a)',' optimal offsets:'
      do i=1,2
         associate(oi=>fdconfig%fdoffsets(i,i,1))
            write(lbl,'(a,i0,a)') 'g(',i,')'
            print'(4x,a16,a,es12.4)',lbl,' ->',oi
         end associate
      end do
      do i=1,2
         do j=1,i
            if(i.eq.j) then
               associate(od=>fdconfig%fdoffsets(i,i,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row/col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',od
               end associate
            else
               associate(&
                     or=>fdconfig%fdoffsets(i,j,1),&
                     oc=>fdconfig%fdoffsets(i,j,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row'
                  print'(4x,a16,a,es12.4)',lbl,' ->',or
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',oc
               end associate
            end if
         end do
      end do
      !================================================================
      !  switch back to central/approx=.false. for grad/Hess demo
      fdconfig%fdscheme='central'
      fdconfig%fdoffsets=&
         get_fdoffsets(&
            func,v,fnoises=fnoises,&
            minoffsets=minoffsets,&
            maxoffsets=maxoffsets,&
            approx=.false.)
      !================================================================
      print'(/,a)',' [3] numerical gradient  (central differences)'
      print'(a)',repeat('-',60)
      print'(a)',' calling get_num_grad ...'
      ncalls=0
      g=get_num_grad(func,v,fixed)
      print'(a,i0)',' function calls   : ',ncalls
      print'(a)',' offsets used:'
      do i=1,2
         associate(oi=>fdconfig%fdoffsets(i,i,1))
            write(lbl,'(a,i0,a)') 'g(',i,')'
            print'(4x,a16,a,es12.4)',lbl,' ->',oi
         end associate
      end do
      print'(a)'
      print'(a,3(2x,a12))',&
            '      ','   numerical','    analytic','    rel. err'
      do i=1,2
         associate(gi=>g(i),ga=>g_ana(i),ge=>g_err(i))
            ge=dabs(gi-ga)/max(dabs(ga),epsilon(0.d0))
            print'(a,i0,a,3(2x,es12.4))',&
                  '  g(',i,')=',gi,ga,ge
         end associate
      end do
      !================================================================
      print'(/,a)',' [4] numerical Hessian  (central differences)'
      print'(a)',repeat('-',60)
      print'(a)',' calling get_num_hess ...'
      ncalls=0
      h=get_num_hess(func,v,fixed)
      print'(a,i0)',' function calls   : ',ncalls
      print'(a)',' offsets used:'
      do i=1,2
         do j=1,i
            if(i.eq.j) then
               associate(od=>fdconfig%fdoffsets(i,i,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row/col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',od
               end associate
            else
               associate(&
                     or=>fdconfig%fdoffsets(i,j,1),&
                     oc=>fdconfig%fdoffsets(i,j,2))
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):row'
                  print'(4x,a16,a,es12.4)',lbl,' ->',or
                  write(lbl,'(a,i0,a,i0,a)') 'H(',i,',',j,'):col'
                  print'(4x,a16,a,es12.4)',lbl,' ->',oc
               end associate
            end if
         end do
      end do
      print'(a)'
      print'(a,3(2x,a12))',&
            '         ','   numerical','    analytic','    rel. err'
      do i=1,2
         do j=1,i
            associate(&
                  hij=>h(i,j),&
                  haj=>h_ana(i,j),&
                  hej=>h_err(i,j))
               if(dabs(haj).lt.epsilon(0.d0)) then
                  print'(a,i0,a,i0,a,2(2x,es12.4),2x,a12)',&
                        '  H(',i,',',j,')=',hij,haj,'         n/a'
               else
                  hej=dabs(hij-haj)/dabs(haj)
                  print'(a,i0,a,i0,a,3(2x,es12.4))',&
                        '  H(',i,',',j,')=',hij,haj,hej
               end if
            end associate
         end do
      end do
# ifdef USE_LAPACK      
      !================================================================
      print'(/,a)',' [5] spectral decomposition of Hessian'
      print'(a)',repeat('-',60)
      print'(a)',' calling diag_num_hess ...'
      call diag_num_hess(func,v,eigvals,eigvects,fixed)
      print'(a)'
      print'(a)',' eigenvalues:'
      do i=1,2
         associate(lami=>eigvals(i))
            print'(4x,a,i0,a,es12.4)','lambda(',i,')=',lami
         end associate
      end do
      print'(a)'
      print'(a)',' eigenvectors (columns):'
      do i=1,2
         associate(vi=>eigvects(:,i))
            print'(4x,a,i0,a,2es12.4)','v',i,'=',vi
         end associate
      end do
      if(eigvals(1).gt.0.d0) then
         print'(/,a)',&
               ' => Hessian is positive definite (local min)'
      else if(eigvals(2).lt.0.d0) then
         print'(/,a)',&
               ' => Hessian is negative definite (local max)'
      else
         print'(/,a)',' => Hessian is indefinite'
      end if
# endif
      print'(a/)',repeat('=',60)
      !***
      contains
      !*****************************************************************
         function func(v) result(res)
      !*****************************************************************
         implicit none
         real(r64) :: &
            v(:),res
         !***
         ncalls=ncalls+1
         res=dsin(a*v(1))*dexp(-b*v(2)**2)
         !***
         end
      !***
      end
