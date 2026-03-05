!***********************************************************************
      program main
!***********************************************************************
      use iso_fortran_env, only : &
         r64=>real64
      use mod_fdbase, only : &
         fdconfig, get_num_grad
      use mod_fdtune, only : &
         get_noise_table, get_fdoffsets
      implicit none
      integer, parameter :: &
         nfun=50
      integer :: &
         ncalls=0
      character(len=5):: &
         f_names(nfun)
      real(r64):: &
         omega,&
         ini_omega,&
         err,&
         g_ana,&
         ratio,&
         v(1),&
         f0,&
         max_err,&
         tmp,&
         eps,&
         p_a,p_b,p_c,p_d,p_e,p_xc,&
         res_list(2*nfun)
      integer:: &
         i,&
         j,&
         mode,&
         imax,&
         n_fd,&
         tmp_idx,&
         res_idx(2*nfun),&
         id,&
         ntot,&
         fd_stats(2*nfun)
      real(r64), allocatable:: &
         noise_table(:,:),&
         fnoises(:),&
         minoffsets(:),&
         maxoffsets(:),&
         g_fd(:)
      !***
      allocate(maxoffsets(1))
      eps=epsilon(0.d0)
      ntot=2*nfun
      do i=1,nfun
         write(f_names(i),'("f_",i2.2)') i
      end do
      max_err=0.d0
      imax=-1
      !---header---
      print '(a)',repeat('=',75)
      print '(a)',&
            '  ID  scheme   func   ratio(r)       g_ana      '//&
            'omega        err   #func'
      print '(a)',repeat('-',75)
      !---main loop---
      do mode=1,2
         if(mode.eq.1) then
            fdconfig%fdscheme='central'
         else
            fdconfig%fdscheme='forward'
         endif
         do i=1,nfun
            ncalls=0
            id=(mode-1)*nfun+i
            call init_params(p_a,p_b,p_c,p_d,p_e,p_xc)
            v(1)=0.5d0
            f0=func(v)
            call set_max_offset(v(1),maxoffsets(1))
            noise_table=get_noise_table(func,v,maxoffsets)
            minoffsets=noise_table(:,1)
            fnoises=noise_table(:,2)
            ncalls=0
            fdconfig%fdoffsets=&
               get_fdoffsets(func,v,fnoises=fnoises,&
                             minoffsets=minoffsets,&
                             maxoffsets=maxoffsets)
            n_fd=ncalls
            omega=fdconfig%fdoffsets(1,1,1)
            ini_omega=&
                dsqrt(dsqrt(minoffsets(1)*maxoffsets(1))*minoffsets(1))
            ratio=omega/ini_omega
            g_fd=get_num_grad(func,v)
            g_ana=fun_deriv(v(1))
            if(dabs(g_ana).lt.eps) then
               err=0.d0
            else
               err=dabs(g_fd(1)-g_ana)/dabs(g_ana)
            endif
            res_list(id)=err
            res_idx(id)=id
            fd_stats(id)=n_fd
            if(err.gt.max_err) then
               max_err=err
               imax=id
            endif
            print '(i4,2x,a7,2x,a5,2x,f8.2,2x,es10.2,2x,es9.1,2x,'//&
                  'es9.1,3x,i5)',&
                  id,fdconfig%fdscheme,f_names(i),ratio,&
                  g_ana,omega,err,n_fd
         end do
         if(mode.eq.1) print '(a)',repeat('-',75)
      end do
      !--- sort for summary ---
      do i=1,ntot-1
         do j=i+1,ntot
            if(res_list(i).gt.res_list(j)) then
               tmp=res_list(i)
               res_list(i)=res_list(j)
               res_list(j)=tmp
               tmp_idx=res_idx(i)
               res_idx(i)=res_idx(j)
               res_idx(j)=tmp_idx
            endif
         end do
      end do
      print '(a)',repeat('=',75)
      print '(" summary (n=",i0,")")',ntot
      print '(a22,es10.2)','median error: ',res_list(ntot/2)
      print '(a22,es10.2,a,i0,a)',&
            'max err: ',max_err,' (ID: ',imax,')'
      print '(a22,f8.2)','avg func calls: ',sum(fd_stats)/dble(ntot)
      print '(a22,i8)','max func calls: ',maxval(fd_stats)
      contains
      !*****************************************************************
         subroutine set_max_offset(x,hmax)
      !*****************************************************************
         real(r64):: &
            x
         real(r64):: &
            hmax
         real(r64):: &
            eps_val
         !***
         eps_val=epsilon(x)
         if(fdconfig%fdscheme.eq.'central') then
            hmax=1.d5*eps_val**(1.d0/3.d0)*max(1.d0,dabs(x))
         else
            hmax=1.d5*dsqrt(eps_val)*max(1.d0,dabs(x))
         endif
         !***
         end
      !*****************************************************************
         subroutine init_params(a,b,c,d,e,xc)
      !*****************************************************************
         implicit none
         real(r64):: &
            a,b,c,d,e,xc,r
         !***
         call random_number(r)
         select case(i)
            case(1)
               a=1.d0+10.d0*r
               b=-5.d0+10.d0*r
            case(2)
               a=0.01d0+0.1d0*r
               b=10.d0+50.d0*r
               c=r
            case(3)
               a=100.d0*r
               b=0.5d0+2.d0*r
            case(4)
               a=1.d-5+1.d-3*r
               b=5.d0+5.d0*r
            case(5)
               a=1.d0+2.d0*r
               b=0.1d0+0.5d0*r
            case(6)
               a=0.5d0*r
               b=1.5d0+3.d0*r
            case(7)
               a=1.d1+1.d2*r
               b=0.1d0+r
               xc=0.73d0
            case(8)
               a=1.d0+r
               b=-2.d0+r
               c=3.d0+r
               d=0.1d0*r
            case(9)
               a=0.5d0+0.5d0*r
               xc=0.73d0+1.d-4*r
            case(10)
               a=2.d0+2.d0*r
               b=10.d0+10.d0*r
            case(11)
               a=1.d0+9.d0*r
               b=0.01d0+0.1d0*r
               xc=0.73d0
            case(12)
               a=r
               b=1.d0+5.d0*r
            case(13)
               a=1.d-2+1.d-1*r
               b=10.d0*r
            case(14)
               a=1.d0+r
               b=0.1d0+0.9d0*r
            case(15)
               a=5.d0+5.d0*r
               b=2.d0+r
            case(16)
               a=0.1d0+0.2d0*r
               b=0.5d0+0.5d0*r
            case(17)
               a=5.d0*r
               b=3.d0*r
            case(18)
               a=0.2d0+0.8d0*r
               xc=0.73d0
            case(19)
               a=1.d2+1.d3*r
               b=1.d0+r
            case(20)
               a=1.d0+r
               xc=0.73d0
            case(21)
               a=1.d-8+1.d-7*r
               b=1.d0+r
            case(22)
               a=1.d8+1.d9*r
               b=r
            case(23)
               a=3.14d0*r
               b=10.d0+90.d0*r
            case(24)
               a=1.d0+r
               b=0.001d0+0.005d0*r
            case(25)
               a=1.d12+1.d13*r
               b=r
            case(26)
               a=1.d-12+1.d-13*r
               b=r
            case(27)
               a=r
               b=1.d0+r
               c=2.d0+r
               d=3.d0+r
               e=4.d0+r
            case(28)
               a=0.001d0+0.002d0*r
               b=50.d0+50.d0*r
            case(29)
               a=10.d0+10.d0*r
               b=1.d-2+1.d-2*r
            case(30)
               a=1.d0
               b=0.5d0
               xc=0.73d0
            case(31)
               a=0.5d0+r
               b=0.1d0+0.2d0*r
            case(32)
               a=2.d0+r
               b=0.05d0+0.05d0*r
            case(33)
               a=10.d0*r
               b=2.d0*r
               c=5.d0+5.d0*r
            case(34)
               a=1.d0+4.d0*r
               b=r
            case(35)
               a=1.5d0+r
               b=10.d0*r
               c=0.7d0+0.1d0*r
            case(36)
               a=1.d0+r
               b=100.d0*r
            case(37)
               a=0.1d0*r
               b=0.5d0+0.5d0*r
            case(38)
               a=1.d0+r
               b=0.2d0+0.3d0*r
            case(39)
               a=0.1d0+r
               b=1.d0+r
               c=2.d0+2.d0*r
            case(40)
               a=10.d0+r
               c=0.05d0
               xc=0.73d0               
            case(41)
               a=1.d0+r
               b=0.1d0+0.1d0*r
            case(42)
               a=2.d0+r
               b=0.05d0+0.05d0*r
            case(43)
               a=0.5d0+r
               b=0.1d0+0.4d0*r
            case(44)
               a=1.d0+r
               b=-1.d0+0.5d0*r
            case(45)
               a=10.d0*r
               b=5.d0+5.d0*r
            case(46)
               a=0.5d0+r
               b=2.d0*r
            case(47)
               a=0.1d0+0.9d0*r
               b=2.d0+8.d0*r
            case(48)
               a=0.5d0*r
               b=0.5d0*r
               xc=0.73d0
            case(49)
               a=1.d0+r
               b=0.1d0+2.d0*r
            case(50)
               a=0.1d0*r
               b=0.2d0*r
               c=0.3d0*r
               d=0.4d0*r
               e=0.5d0*r
               xc=0.73d0
            case default
               a=1.d0
               b=1.d0
         end select
         !***
         end
      !*****************************************************************
         function func(v) result(res)
      !*****************************************************************
         implicit none
         real(r64):: &
            v(:),res,x
         !***
         ncalls=ncalls+1
         x=v(1)
         select case(i)
            case(1)
               res=p_a*x+p_b
            case(2)
               res=p_a*x**2+p_b*x+p_c
            case(3)
               res=p_a*dsin(p_b*x)
            case(4)
               res=p_a*dexp(p_b*x)
            case(5)
               res=p_a*dlog(x+p_b)
            case(6)
               res=p_a/(x+p_b)
            case(7)
               res=p_a*dtanh(p_b*(x-p_xc))
            case(8)
               res=p_a*x**3+p_b*x**2+p_c*x+p_d
            case(9)
               res=p_a*dsqrt(dabs(x-p_xc)+1.d-14)
            case(10)
               res=p_a*dcos(p_b*x**2)
            case(11)
               res=p_a*dexp(-p_b/(max(1.d-12,(x-p_xc)**2)+0.1d0))
            case(12)
               res=p_a/(1.d0+dexp(-p_b*x))
            case(13)
               res=p_a*x**4+p_b*x**2
            case(14)
               res=p_a*dsin(p_b*x)/max(1.d-12,x)
            case(15)
               res=p_a*dlog10(x+p_b)
            case(16)
               res=p_a*x*dexp(p_b*x)
            case(17)
               res=p_a*dsin(x)+p_b*dcos(x)
            case(18)
               res=p_a*dabs(x-(p_xc+0.5d0))**(1.d0/3.d0)
            case(19)
               res=p_a/(x+p_b)**2
            case(20)
               res=p_a*dabs(x-p_xc)**1.5d0
            case(21)
               res=p_a*x**2+p_b*x+1.d-5
            case(22)
               res=p_a*x**2-p_b*x
            case(23)
               res=p_a*dcos(p_b*x)
            case(24)
               res=p_a*dtan(p_b*x)
            case(25)
               res=p_a*dexp(-x)
            case(26)
               res=p_a*dlog(x**2+1.d0)
            case(27)
               res=p_a*x**4+p_b*x**3+p_c*x**2+p_d*x+p_e
            case(28)
               res=p_a*dexp(p_b*x)
            case(29)
               res=p_a*dsin(x/max(1.d-12,p_b))   
            case(30)
               res=p_a*dexp(-p_b*(x-p_xc)**2)*dsin(x)
            case(31)
               res=p_a*dsinh(p_b*x)
            case(32)
               res=p_a*dcosh(p_b*x)
            case(33)
               res=p_a*x+p_b*dsin(p_c*x)
            case(34)
               res=p_a/(1.d0+x**2)+p_b*x
            case(35)
               res=p_a*datan(p_b*(x-p_c))
            case(36)
               res=p_a*dexp(-p_b*x**2)
            case(37)
               res=p_a*x**5+p_b*x**3
            case(38)
               res=p_a*dtan(p_b*x)
            case(39)
               res=p_a*x**p_c+p_b
            case(40)
               res=p_a*x**2*dtanh((x-p_xc)/p_c)
            case(41)
               res=p_a*dacos(0.1d0*p_b*x)
            case(42)
               res=p_a*dasin(0.1d0*p_b*x)
            case(43)
               res=p_a*dexp(p_b*dsqrt(x))
            case(44)
               res=p_a*dcos(x)*dexp(p_b*x)
            case(45)
               res=p_a/(x**2+p_b*x+1.d0)
            case(46)
               res=p_a*x+p_b*dsin(x)**2
            case(47)
               res=p_a*dlog(dabs(dcos(p_b*x))+0.1d0)
            case(48)
               res=p_a*(x-p_xc)**2*dtanh((x-p_xc)/0.05d0)               
            case(49)
               res=p_a*x**(-0.5d0)+p_b*x
            case(50)
               res=p_a*x**4+p_b*x**3+p_c*x**2+p_d*x+p_e
            case default
               res=p_a*x
         end select
         !***
         end
      !*****************************************************************
         function fun_deriv(x) result(res)
      !*****************************************************************
         implicit none
         real(r64):: &
            x,res,tmp,s
         !***
         select case(i)
            case(1)
               res=p_a
            case(2)
               res=2.d0*p_a*x+p_b
            case(3)
               res=p_a*p_b*dcos(p_b*x)
            case(4)
               res=p_a*p_b*dexp(p_b*x)
            case(5)
               res=p_a/(x+p_b)
            case(6)
               res=-p_a/(x+p_b)**2
            case(7)
               res=p_a*p_b*(1.d0-dtanh(p_b*(x-p_xc))**2)
            case(8)
               res=3.d0*p_a*x**2+2.d0*p_b*x+p_c
            case(9)
               s=merge(1.d0,-1.d0,x.gt.p_xc)
               res=0.5d0*p_a*s/dsqrt(dabs(x-p_xc)+1.d-14)
            case(10)
               res=-2.d0*p_a*p_b*x*dsin(p_b*x**2)
            case(11)
               tmp=(x-p_xc)**2+0.1d0
               res=func((/x/))*(2.d0*p_b*(x-p_xc)/tmp**2)
            case(12)
               tmp=dexp(-p_b*x)
               res=p_a*p_b*tmp/(1.d0+tmp)**2
            case(13)
               res=4.d0*p_a*x**3+2.d0*p_b*x
            case(14)
               res=p_a*(p_b*x*dcos(p_b*x)-dsin(p_b*x))/x**2
            case(15)
               res=p_a/((x+p_b)*dlog(10.d0))
            case(16)
               res=p_a*dexp(p_b*x)*(1.d0+p_b*x)
            case(17)
               res=p_a*dcos(x)-p_b*dsin(x)
            case(18)
               tmp=x-(p_xc+0.5d0)
               s=dsign(1.d0,tmp)
               res=(p_a/3.d0)*dabs(tmp)**(-2.d0/3.d0)*s               
            case(19)
               res=-2.d0*p_a/(x+p_b)**3
            case(20)
               s=merge(1.d0,-1.d0,x.gt.p_xc)
               res=1.5d0*p_a*dsqrt(dabs(x-p_xc))*s
            case(21)
               res=2.d0*p_a*x+p_b
            case(22)
               res=2.d0*p_a*x-p_b
            case(23)
               res=-p_a*p_b*dsin(p_b*x)
            case(24)
               res=p_a*p_b/(dcos(p_b*x)**2)
            case(25)
               res=-p_a*dexp(-x)
            case(26)
               res=2.d0*p_a*x/(x**2+1.d0)
            case(27)
               res=4.d0*p_a*x**3+3.d0*p_b*x**2+2.d0*p_c*x+p_d
            case(28)
               res=p_a*p_b*dexp(p_b*x)
            case(29)
               res=(p_a/p_b)*dcos(x/p_b)               
            case(30)
               res=p_a*dexp(-p_b*(x-p_xc)**2)*(dcos(x)-&
                  2.d0*p_b*(x-p_xc)*dsin(x))
            case(31)
               res=p_a*p_b*dcosh(p_b*x)
            case(32)
               res=p_a*p_b*dsinh(p_b*x)
            case(33)
               res=p_a+p_b*p_c*dcos(p_c*x)
            case(34)
               res=-2.d0*p_a*x/(1.d0+x**2)**2+p_b
            case(35)
               res=p_a*p_b/(1.d0+(p_b*(x-p_c))**2)
            case(36)
               res=-2.d0*p_a*p_b*x*dexp(-p_b*x**2)
            case(37)
               res=5.d0*p_a*x**4+3.d0*p_b*x**2
            case(38)
               res=p_a*p_b/(dcos(p_b*x)**2)
            case(39)
               res=p_a*p_c*x**(p_c-1.d0)
            case(40)
               res=2.d0*p_a*x*dtanh((x-p_xc)/p_c)+&
                  p_a*x**2*(1.d0-dtanh((x-p_xc)/p_c)**2)/p_c
            case(41)
               res=-0.1d0*p_a*p_b/dsqrt(1.d0-(0.1d0*p_b*x)**2)
            case(42)
               res=0.1d0*p_a*p_b/dsqrt(1.d0-(0.1d0*p_b*x)**2)
            case(43)
               res=p_a*p_b*dexp(p_b*dsqrt(x))/(2.d0*dsqrt(x))
            case(44)
               res=p_a*dexp(p_b*x)*(p_b*dcos(x)-dsin(x))
            case(45)
               res=-p_a*(2.d0*x+p_b)/(x**2+p_b*x+1.d0)**2
            case(46)
               res=p_a+2.d0*p_b*dsin(x)*dcos(x)
            case(47)
               res=-p_a*p_b*dtan(p_b*x)/(1.d0+0.1d0/dabs(dcos(p_b*x)))
            case(48)
               res=p_a*(2.d0*(x-p_xc)*dtanh((x-p_xc)/0.05d0)+&
                  (x-p_xc)**2*(1.d0-dtanh((x-p_xc)/0.05d0)**2)/0.05d0)
            case(49)
               res=-0.5d0*p_a*x**(-1.5d0)+p_b
            case(50)
               res=4.d0*p_a*x**3+3.d0*p_b*x**2+2.d0*p_c*x+p_d
            case default
               res=p_a
         end select
         !***
         end
      end
