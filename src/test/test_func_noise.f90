!***********************************************************************
      program main
!***********************************************************************
      use iso_fortran_env, only : &
         r64=>real64
      use mod_fdtune, only : &
         get_noise_table
      implicit none
      integer, parameter :: &
         nfun=8,&
         ntot=150
      integer :: &
         ncalls=0
      real(r64) :: &
         v(1),&
         f0,&
         t_nu,&
         est_nu,&
         fac,&
         tmp,&
         abs_max,&
         r,&
         p_a,p_w,p_nl,p_xc,&
         res_list(ntot)
      integer :: &
         i,j,&
         imax,&
         test_idx,&
         stats(ntot)
      real(r64), allocatable :: &
         noise_table(:,:),&
         fnoises(:),&
         maxoffsets(:)
      character(len=8), parameter :: &
         f_names(nfun)=(/&
            'g_osc   ',&
            'g_peak  ',&
            'g_gauss ',&
            'g_corner',&
            'g_cont  ',&
            'g_disc  ',&
            'rosen   ',&
            'sing    '/)
      !***
      maxoffsets=[1.d0]
      abs_max=-1.d0
      imax=-1
      !--- header ---
      print'(a)',repeat('=',58)
      print'(a4,2x,a,6x,a10,2x,a10,2x,a8,2x,a5)',&
            ' ID','func','nu(true)','nu(est)','ratio','#func'
      print'(a)',repeat('-',58)
      !--- main loop ---
      do i=1,ntot
         test_idx=mod(i-1,nfun)
         call random_number(r)
         p_nl=1.d1**(-15.d0+13.d0*r)
         if(test_idx.ge.6) then
            call random_number(r)
            v(1)=-0.5d0+2.d0*r
         else
            call random_number(r)
            p_a=0.1d0+99.9d0*r
            call random_number(r)
            p_w=r
            call random_number(r)
            p_xc=r
            call random_number(r)
            v(1)=0.1d0+0.8d0*r
         endif
         ncalls=0
         noise_table=get_noise_table(func,v,maxoffsets)
         fnoises=noise_table(:,2)
         est_nu=fnoises(1)
         f0=f_exact(test_idx,v(1))
         t_nu=max(p_nl,epsilon(0.d0)*max(1.d0,dabs(f0)))
         fac=max(est_nu/t_nu,t_nu/est_nu)
         res_list(i)=fac
         if(fac.gt.abs_max) then
            abs_max=fac
            imax=i
         endif
         stats(i)=ncalls
         print'(i4,2x,a8,2x,es10.2,2x,es10.2,2x,f8.2,2x,i5)',&
               i,f_names(test_idx+1),t_nu,est_nu,fac,ncalls
      end do
      !--- sort for summary ---
      do i=1,ntot-1
         do j=i+1,ntot
            if(res_list(i).gt.res_list(j)) then
               tmp=res_list(i)
               res_list(i)=res_list(j)
               res_list(j)=tmp
            endif
         end do
      end do
      print'(a)',repeat('=',58)
      print'(" summary (n=",i0,")")',ntot
      print'(a22,f8.2)','median fac: ',res_list(ntot/2)
      print'(a22,f8.2,a,i0,a)','max fac: ',abs_max,' (ID: ',imax,')'
      print'(a22,f8.1)','avg func calls: ',sum(stats)/dble(ntot)
      print'(a22,i8)','max func calls: ',maxval(stats)
      !***
      contains
      !*****************************************************************
         function rand_normal() result(res)
      !*****************************************************************
         implicit none
         real(r64) :: &
            res,&
            u(2)
         !***
         call random_number(u)
         res=dsqrt(-2.d0*dlog(u(1)))*dcos(8.d0*datan(1.d0)*u(2))
         !***
         end
      !*****************************************************************
         function f_exact(idx,x) result(res)
      !*****************************************************************
         implicit none
         integer :: idx
         real(r64) :: &
            x,&
            res
         !***
         select case(idx)
            case(0)
               res=dcos(2.d0*3.1415926535d0*p_w+p_a*x)
            case(1)
               res=1.d0/(p_a**(-2.d0)+(x-p_w)**2)
            case(2)
               res=dexp(-(p_a**2*(x-p_w)**2))
            case(3)
               res=dexp(-p_a*dabs(x-p_w))
            case(4)
               if(x.ge.p_w) then
                  res=dexp(-p_a*(x-p_w))
               else
                  res=dexp(p_a*(x-p_w))
               endif
            case(5)
               if(x.lt.p_w) then
                  res=0.d0
               else
                  res=dexp(p_a*x)
               endif
            case(6)
               res=(1.d0-x)**2+100.d0*(x**2-x)**2
            case(7)
               res=1.d0/(dabs(x-0.1d0)**0.5d0+0.01d0)
         end select
         !***
         end
      !*****************************************************************
         function func(v) result(res)
      !*****************************************************************
         implicit none
         real(r64) :: &
            v(:),&
            res
         !***
         ncalls=ncalls+1
         res=f_exact(test_idx,v(1))+p_nl*rand_normal()
         !***
         end
      !***
      end
