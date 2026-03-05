!***********************************************************************
      program main
!***********************************************************************
      use iso_fortran_env, only : &
         r64=>real64
      use mod_fdbase, only : &
         fdconfig, get_num_grad, get_num_hess
      use mod_fdtune, only : &
         get_noise_table, get_fdoffsets
      implicit none
      logical :: approx
      integer, parameter :: &
         nfun=50,&
         inds(3,2)=reshape((/1,2,2, 1,1,2/),(/3,2/))
      integer :: &
         ncalls=0
      character(len=5) :: &
         f_names(nfun)
      character(len=16) :: &
         arg
      real(r64), allocatable :: &
         noise_table(:,:),&
         fnoises(:),&
         minoffsets(:),&
         maxoffsets(:),&
         hess_fd(:,:),&
         g(:),&
         res_list(:)
      real(r64) :: &
         omega(3),&
         err(3),&
         h_ana(3),&
         ratio(3),&
         v(2),&
         f0,&
         hnj,&
         max_err,&
         tmp,&
         eps,&
         omega0,&
         p_a,p_b,p_c
      integer :: &
         i,j,&
         id,&
         imax,&
         mode,&
         n_fd,&
         ntot,&
         tmp_idx
      integer, allocatable :: &
         res_idx(:),&
         fd_stats(:)
      !***
      !***
      approx=.false.
      if(iargc().gt.0) then
         call getarg(1,arg)
         if(trim(arg).eq.'approx') approx=.true.
      endif
      allocate(maxoffsets(2))
      eps=epsilon(0.d0)
      ntot=merge(1,2,approx)*nfun
      allocate(res_list(ntot),res_idx(ntot),fd_stats(ntot))
      do i=1,nfun
         write(f_names(i),'("f_",i2.2)') i
      end do
      max_err=-1.d0
      imax=-1
      !---header---
      print'(a)',repeat('=',159)
      print'(a)','  ID  scheme   func'//&
           '     r(1,1)    r(2,1)    r(2,2)'//&
           '  h_ana(1,1)  h_ana(2,1)  h_ana(2,2)'//&
           ' omega(1,1) omega(2,1) omega(2,2)'//&
           '   err(1,1)   err(2,1)   err(2,2)  #func'
      print'(a)',repeat('-',159)
      !---main loop---
      do mode=1,2
         if(mode.eq.2.and.approx) exit 
         if(mode.eq.1) then
            fdconfig%fdscheme='central'
         else
            fdconfig%fdscheme='forward'
         endif
         do i=1,nfun
            ncalls=0
            id=(mode-1)*nfun+i
            call init_params(p_a,p_b,p_c)
            v(1)=0.5d0
            v(2)=0.5d0
            f0=func(v)
            do j=1,2
               call set_max_offset(v(j),maxoffsets(j))
            end do
            noise_table=get_noise_table(func,v,maxoffsets)
            minoffsets=noise_table(:,1)
            fnoises=noise_table(:,2)
            ncalls=0
            fdconfig%fdoffsets=&
                  get_fdoffsets(&
                     func,v,fnoises=fnoises,&
                     minoffsets=minoffsets,&
                     maxoffsets=maxoffsets,approx=approx)
            n_fd=ncalls
            hess_fd=get_num_hess(func,v)
            g=get_num_grad(func,v)
            !--- collect per-element results ---
            omega0=dsqrt(dsqrt(minoffsets(1)*maxoffsets(1))*&
                               minoffsets(1))
            do j=1,3
               associate(&
                     ii     => inds(j,1),&
                     jj     => inds(j,2),&
                     haj    => h_ana(j),&
                     errj   => err(j),&
                     omegaj => omega(j),&
                     ratioj => ratio(j))
                  hnj=hess_fd(ii,jj)
                  omegaj=fdconfig%fdoffsets(ii,jj,2)
                  haj=fun_hess(v(1),v(2),ii,jj)
                  if(dabs(haj).lt.eps) then
                     errj=0.d0
                  else
                     errj=dabs(haj-hnj)/dabs(haj)
                  endif
                  ratioj=omegaj/omega0
               end associate
            end do
            !--- track worst error across all elements ---
            do j=1,3
               if(err(j).gt.max_err) then
                  max_err=err(j)
                  imax=id
               endif
            end do
            res_list(id)=maxval(err)
            res_idx(id)=id
            fd_stats(id)=n_fd
            print'(i4,2x,a7,2x,a5,&
                   &3(2x,f8.2),3(2x,es10.2),&
                   &3(2x,es9.1),3(2x,es9.1),2x,i5)',&
                  id,fdconfig%fdscheme,f_names(i),&
                  ratio(1),ratio(2),ratio(3),&
                  h_ana(1),h_ana(2),h_ana(3),&
                  omega(1),omega(2),omega(3),&
                  err(1),err(2),err(3),n_fd
         end do
         if(mode.eq.1) print'(a)',repeat('-',159)
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
      print'(a)',repeat('=',159)
      print'(" summary (n=",i0,")")',ntot
      print'(a22,es10.2)','median error: ',res_list(ntot/2)
      print'(a22,es10.2,a,i0,a)',&
            'max err: ',max_err,' (ID: ',imax,')'
      print'(a22,f8.2)','avg func calls: ',&
            sum(fd_stats)/dble(ntot)
      print'(a22,i8)','max func calls: ',maxval(fd_stats)
      contains
      !*****************************************************************
         subroutine set_max_offset(x,hmax)
      !*****************************************************************
         implicit none
         real(r64) :: &
            x,hmax,eps_val
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
         subroutine init_params(a,b,c)
      !*****************************************************************
         implicit none
         real(r64) :: &
            a,b,c,r
         !***
         call random_number(r)
         select case(i)
            case(1)
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(2)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(3)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(4)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(5)        
               a=1.d0+9.d0*r
            case(6)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(7)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(8)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(9)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(10)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(11)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(12)
               a=1.d0+9.d0*r
               b=0.5d0+0.5d0*r
            case(13)
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(14)
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(15)        
               a=1.d0+9.d0*r
            case(16)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(17)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(18)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(19)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(20)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(21)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(22)        
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(23)        
               a=1.d0+9.d0*r
            case(24)        
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(25)
               a=1.d0+9.d0*r
               b=0.5d0+0.5d0*r
            case(26)
               a=1.d0+9.d0*r
            case(27)
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(28)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(29)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(30)       
               a=1.d0+9.d0*r
            case(31)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+9.d0*r
            case(32)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(33)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+9.d0*r
            case(34)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(35)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(36)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(37)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(38)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(39)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(40)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+4.d0*r
            case(41)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+4.d0*r
            case(42)       
               a=1.d0+9.d0*r
            case(43)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(44)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(45)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
            case(46)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
               c=1.d0+9.d0*r
            case(47)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case(48)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(49)       
               a=1.d0+9.d0*r
               b=1.d0+9.d0*r
            case(50)       
               a=1.d0+9.d0*r
               b=1.d0+4.d0*r
               c=1.d0+4.d0*r
            case default
               a=1.d0
               b=1.d0
               c=1.d0
         end select
         !***
         end
      !*****************************************************************
         function func(v) result(res)
      !*****************************************************************
         real(r64) :: v(:),res,x,y
         !***
         ncalls=ncalls+1
         x=v(1)
         y=v(2)
         select case(i)
            case(1)
               res=p_a*dexp(p_b*x)*dsin(p_c*y)
            case(2)
               res=p_a*x**2+p_b*y**4
            case(3)
               res=p_a*dsin(p_b*x)*dcos(p_c*y)
            case(4)
               res=p_a*dtanh(p_b*x)*dexp(p_c*y)
            case(5)
               res=p_a*x**3*y**2
            case(6)
               res=p_a*dexp(p_b*x**2)*dsin(p_c*y)
            case(7)
               res=p_a*dlog(x**2+1.d0)+p_b*y**4
            case(8)
               res=p_a*x**2*y+p_b*y**3
            case(9)
               res=p_a/(x**2+1.d0)*dexp(p_b*y)
            case(10)
               res=p_a*dsin(p_b*x**2)*dcos(p_c*y)
            case(11)
               res=p_a*dcos(p_b*x)+p_c*y**3
            case(12)
               res=p_a*dexp(p_b*x*y)
            case(13)
               res=p_a*dsin(p_b*x+p_c*y)
            case(14)
               res=p_a*x**4+p_b*y**2
            case(15)
               res=p_a/(x**2+y**2+1.d0)
            case(16)
               res=p_a*dlog(p_b*x**2+p_c*y**2+1.d0)
            case(17)
               res=p_a*x**2*dexp(p_b*y)
            case(18)
               res=p_a*dtanh(p_b*x)*dtanh(p_c*y)
            case(19)
               res=p_a*dcos(p_b*x**2+p_c*y)
            case(20)
               res=p_a*x**3+p_b*y**3
            case(21)
               res=p_a*dexp(p_b*x)*dcos(p_c*y**2)
            case(22)
               res=p_a*dsinh(p_b*x)*dcosh(p_c*y)
            case(23)
               res=p_a*x**2/(y**2+1.d0)
            case(24)
               res=p_a*dsin(p_b*x)*dsinh(p_c*y)
            case(25)
               res=p_a*dexp(p_b*(x**2+y**2))
            case(26)
               res=p_a*x**2*y**3
            case(27)
               res=p_a*dcos(p_b*x)*dsin(p_c*y**2)
            case(28)
               res=p_a/(p_b*x+1.d0)*dexp(p_c*y)
            case(29)
               res=p_a*dlog(p_b*x**2+1.d0)*y**2
            case(30)
               res=p_a*x**4*y**2
            case(31)
               res=p_a*dtanh(p_b*x**2)*dsin(p_c*y)
            case(32)
               res=p_a*dexp(p_b*x)*y**4
            case(33)
               res=p_a*dsin(p_b*x**3)*dcos(p_c*y)
            case(34)
               res=p_a*dlog(x**2+1.d0)*dsin(p_b*y)
            case(35)
               res=p_a*x**2*dcos(p_b*y**2)
            case(36)
               res=p_a*dexp(p_b*x**2)*dcosh(p_c*y)
            case(37)
               res=p_a*dsin(p_b*x)/(y**2+1.d0)
            case(38)
               res=p_a*x**3*dexp(p_b*y)
            case(39)
               res=p_a*dtanh(p_b*x)*y**4
            case(40)
               res=p_a*dcos(p_b*x)*dexp(p_c*y**2)
            case(41)
               res=p_a*dsin(p_b*x**2)*dsinh(p_c*y)
            case(42)
               res=p_a*x**3*dlog(y**2+1.d0)
            case(43)
               res=p_a*dexp(p_b*x)*dtanh(p_c*y**2)
            case(44)
               res=p_a*dcos(p_b*x**3)*y**2
            case(45)
               res=p_a*x**2*dsinh(p_b*y)
            case(46)
               res=p_a*dlog(p_b*x**2+1.d0)*dcos(p_c*y)
            case(47)
               res=p_a*dexp(p_b*x**2)*dtanh(p_c*y)
            case(48)
               res=p_a*dsin(p_b*x)*y**3
            case(49)
               res=p_a*x**4*dcos(p_b*y)
            case(50)
               res=p_a*dtanh(p_b*x**3)*dexp(p_c*y)
            case default
               res=p_a*x*y
         end select
         !***
         end
      !*****************************************************************
         function fun_hess(x,y,ii,jj) result(res)
      !*****************************************************************
         implicit none
         integer :: &
            ii,jj
         real(r64) :: &
            x,y,res,tx,ty,&
            sx,sy,ex,ey,t1,t2
         !***
         select case(i)
            case(1)
               ex=dexp(p_b*x)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*p_b**2*ex*dsin(p_c*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*p_c*ex*dcos(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*ex*dsin(p_c*y)
            case(2)
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a
               if(ii.eq.2.and.jj.eq.1) res=0.d0
               if(ii.eq.2.and.jj.eq.2) res=12.d0*p_b*y**2
            case(3)
               sx=dsin(p_b*x)
               sy=dcos(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*sx*sy
               if(ii.eq.2.and.jj.eq.1) &
                  res=-p_a*p_b*p_c*dcos(p_b*x)*dsin(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*sx*sy
            case(4)
               tx=dtanh(p_b*x)
               ey=dexp(p_c*y)
               t1=1.d0-tx**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b**2*tx*t1*ey
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*p_c*t1*ey
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*tx*ey
            case(5)
               if(ii.eq.1.and.jj.eq.1) res=6.d0*p_a*x*y**2
               if(ii.eq.2.and.jj.eq.1) res=6.d0*p_a*x**2*y
               if(ii.eq.2.and.jj.eq.2) res=2.d0*p_a*x**3
            case(6)
               ex=dexp(p_b*x**2)
               sy=dsin(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*ex*sy*(2.d0*p_b+4.d0*p_b**2*x**2)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*x*ex*dcos(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*ex*sy
            case(7)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*(2.d0-2.d0*x**2)/(x**2+1.d0)**2
               if(ii.eq.2.and.jj.eq.1) res=0.d0
               if(ii.eq.2.and.jj.eq.2) res=12.d0*p_b*y**2
            case(8)
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a*y
               if(ii.eq.2.and.jj.eq.1) res=2.d0*p_a*x
               if(ii.eq.2.and.jj.eq.2) res=6.d0*p_b*y
            case(9)
               ex=dexp(p_b*y)
               t1=x**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*ex*(6.d0*x**2-2.d0)/t1**3
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*x*ex/t1**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_b**2*ex/t1
            case(10)
               t1=p_b*x**2
               sx=dsin(t1)
               sy=dcos(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*sy*&
                  (2.d0*p_b*dcos(t1)-4.d0*p_b**2*x**2*sx)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*x*dcos(t1)*dsin(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*sx*sy
            case(11)
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dcos(p_b*x)
               if(ii.eq.2.and.jj.eq.1) res=0.d0
               if(ii.eq.2.and.jj.eq.2) res=6.d0*p_c*y
            case(12)
               ex=dexp(p_b*x*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*p_b**2*y**2*ex
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*ex*(1.d0+p_b*x*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_b**2*x**2*ex
            case(13)
               t1=p_b*x+p_c*y
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dsin(t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-p_a*p_b*p_c*dsin(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*dsin(t1)
            case(14)
               if(ii.eq.1.and.jj.eq.1) res=12.d0*p_a*x**2
               if(ii.eq.2.and.jj.eq.1) res=0.d0
               if(ii.eq.2.and.jj.eq.2) res=2.d0*p_b
            case(15)
               t1=x**2+y**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=2.d0*p_a*(4.d0*x**2/t1-1.d0)/t1**2
               if(ii.eq.2.and.jj.eq.1) &
                  res=8.d0*p_a*x*y/t1**3
               if(ii.eq.2.and.jj.eq.2) &
                  res=2.d0*p_a*(4.d0*y**2/t1-1.d0)/t1**2
            case(16)
               t1=p_b*x**2+p_c*y**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*(2.d0*p_b*t1-4.d0*p_b**2*x**2)/t1**2
               if(ii.eq.2.and.jj.eq.1) &
                  res=-4.d0*p_a*p_b*p_c*x*y/t1**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*(2.d0*p_c*t1-4.d0*p_c**2*y**2)/t1**2
            case(17)
               ey=dexp(p_b*y)
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a*ey
               if(ii.eq.2.and.jj.eq.1) res=2.d0*p_a*p_b*x*ey
               if(ii.eq.2.and.jj.eq.2) res=p_a*p_b**2*x**2*ey
            case(18)
               tx=dtanh(p_b*x)
               ty=dtanh(p_c*y)
               t1=1.d0-tx**2
               t2=1.d0-ty**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b**2*tx*t1*ty
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*p_c*t1*t2
               if(ii.eq.2.and.jj.eq.2) &
                  res=-2.d0*p_a*p_c**2*tx*ty*t2
            case(19)
               t1=p_b*x**2+p_c*y
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*&
                  (2.d0*p_b*x**2*dcos(t1)+dsin(t1))
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*x*dcos(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*dcos(t1)
            case(20)
               if(ii.eq.1.and.jj.eq.1) res=6.d0*p_a*x
               if(ii.eq.2.and.jj.eq.1) res=0.d0
               if(ii.eq.2.and.jj.eq.2) res=6.d0*p_b*y
            case(21)
               ex=dexp(p_b*x)
               t1=p_c*y**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*p_b**2*ex*dcos(t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*y*ex*dsin(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-2.d0*p_a*p_c*ex*&
                  (2.d0*p_c*y**2*dcos(t1)+dsin(t1))
            case(22)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*p_b**2*dsinh(p_b*x)*dcosh(p_c*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*p_c*dcosh(p_b*x)*dsinh(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*dsinh(p_b*x)*dcosh(p_c*y)
            case(23)
               t1=y**2+1.d0
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a/t1
               if(ii.eq.2.and.jj.eq.1) res=-4.d0*p_a*x*y/t1**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*x**2*(6.d0*y**2-2.d0)/t1**3
            case(24)
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dsin(p_b*x)*dsinh(p_c*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=p_a*p_b*p_c*dcos(p_b*x)*dcosh(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*dsin(p_b*x)*dsinh(p_c*y)
            case(25)
               ex=dexp(p_b*(x**2+y**2))
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*ex*(2.d0*p_b+4.d0*p_b**2*x**2)
               if(ii.eq.2.and.jj.eq.1) &
                  res=4.d0*p_a*p_b**2*x*y*ex
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*ex*(2.d0*p_b+4.d0*p_b**2*y**2)
            case(26)
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a*y**3
               if(ii.eq.2.and.jj.eq.1) res=6.d0*p_a*x*y**2
               if(ii.eq.2.and.jj.eq.2) res=6.d0*p_a*x**2*y
            case(27)
               t1=p_c*y**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dcos(p_b*x)*dsin(t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*y*dsin(p_b*x)*dcos(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*dcos(p_b*x)*&
                     (2.d0*p_c*dcos(t1)-4.d0*p_c**2*y**2*dsin(t1))
            case(28)
               t1=p_b*x+1.d0
               ey=dexp(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=2.d0*p_a*p_b**2*ey/t1**3
               if(ii.eq.2.and.jj.eq.1) &
                  res=-p_a*p_b*p_c*ey/t1**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*ey/t1
            case(29)
               t1=p_b*x**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*y**2*(2.d0*p_b*t1-&
                     4.d0*p_b**2*x**2)/t1**2
               if(ii.eq.2.and.jj.eq.1) &
                  res=4.d0*p_a*p_b*x*y/t1
               if(ii.eq.2.and.jj.eq.2) &
                  res=2.d0*p_a*dlog(t1)
            case(30)
               if(ii.eq.1.and.jj.eq.1) res=12.d0*p_a*x**2*y**2
               if(ii.eq.2.and.jj.eq.1) res=8.d0*p_a*x**3*y
               if(ii.eq.2.and.jj.eq.2) res=2.d0*p_a*x**4
            case(31)
               tx=dtanh(p_b*x**2)
               t1=1.d0-tx**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*dsin(p_c*y)*&
                     (2.d0*p_b*t1-8.d0*p_b**2*x**2*tx*t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*x*t1*dcos(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*tx*dsin(p_c*y)
            case(32)
               ex=dexp(p_b*x)
               if(ii.eq.1.and.jj.eq.1) res=p_a*p_b**2*ex*y**4
               if(ii.eq.2.and.jj.eq.1) res=4.d0*p_a*p_b*ex*y**3
               if(ii.eq.2.and.jj.eq.2) res=12.d0*p_a*ex*y**2
            case(33)
               t1=p_b*x**3
               sy=dcos(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*sy*(6.d0*p_b*x*dcos(t1)-&
                     9.d0*p_b**2*x**4*dsin(t1))
               if(ii.eq.2.and.jj.eq.1) &
                  res=-3.d0*p_a*p_b*p_c*x**2*dcos(t1)*dsin(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*dsin(t1)*sy
            case(34)
               t1=x**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*(2.d0*x**2/t1-1.d0)*dsin(p_b*y)/t1
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*x*dcos(p_b*y)/t1
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_b**2*dlog(t1)*dsin(p_b*y)
            case(35)
               t1=p_b*y**2
               if(ii.eq.1.and.jj.eq.1) res=2.d0*p_a*dcos(t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-4.d0*p_a*p_b*x*y*dsin(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-2.d0*p_a*p_b*x**2*&
                  (2.d0*p_b*y**2*dcos(t1)+dsin(t1))
            case(36)
               ex=dexp(p_b*x**2)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*dcosh(p_c*y)*ex*(2.d0*p_b+4.d0*p_b**2*x**2)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*x*ex*dsinh(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*ex*dcosh(p_c*y)
            case(37)
               t1=y**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dsin(p_b*x)/t1
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*y*dcos(p_b*x)/t1**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*dsin(p_b*x)*(6.d0*y**2-2.d0)/t1**3
            case(38)
               ey=dexp(p_b*y)
               if(ii.eq.1.and.jj.eq.1) res=6.d0*p_a*x*ey
               if(ii.eq.2.and.jj.eq.1) res=3.d0*p_a*p_b*x**2*ey
               if(ii.eq.2.and.jj.eq.2) res=p_a*p_b**2*x**3*ey
            case(39)
               tx=dtanh(p_b*x)
               t1=1.d0-tx**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b**2*tx*t1*y**4
               if(ii.eq.2.and.jj.eq.1) &
                  res=4.d0*p_a*p_b*t1*y**3
               if(ii.eq.2.and.jj.eq.2) res=12.d0*p_a*tx*y**2
            case(40)
               ex=dexp(p_c*y**2)
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dcos(p_b*x)*ex
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*y*dsin(p_b*x)*ex
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*dcos(p_b*x)*ex*(2.d0*p_c+4.d0*p_c**2*y**2)
            case(41)
               t1=p_b*x**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*(2.d0*p_b*x**2*dsin(t1)-&
                  dcos(t1))*dsinh(p_c*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*x*dcos(t1)*dcosh(p_c*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*dsin(t1)*dsinh(p_c*y)
            case(42)
               t1=y**2+1.d0
               if(ii.eq.1.and.jj.eq.1) res=6.d0*p_a*x*dlog(t1)
               if(ii.eq.2.and.jj.eq.1) &
                  res=6.d0*p_a*x**2*y/t1
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*x**3*(2.d0*t1-4.d0*y**2)/t1**2
            case(43)
               ex=dexp(p_b*x)
               t1=p_c*y**2
               ty=dtanh(t1)
               t2=1.d0-ty**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*p_b**2*ex*ty
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*y*ex*t2
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*ex*(2.d0*p_c*t2-8.d0*p_c**2*y**2*ty*t2)
            case(44)
               t1=p_b*x**3
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*y**2*(6.d0*p_b*x*(-dsin(t1))+&
                     (-dcos(t1))*9.d0*p_b**2*x**4)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*3.d0*p_b*x**2*y*dsin(t1)
               if(ii.eq.2.and.jj.eq.2) &
                  res=2.d0*p_a*dcos(t1)
            case(45)
               if(ii.eq.1.and.jj.eq.1) &
                  res=2.d0*p_a*dsinh(p_b*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*x*dcosh(p_b*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_b**2*x**2*dsinh(p_b*y)
            case(46)
               t1=p_b*x**2+1.d0
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*dcos(p_c*y)*(2.d0*p_b*t1-&
                     4.d0*p_b**2*x**2)/t1**2
               if(ii.eq.2.and.jj.eq.1) &
                  res=-2.d0*p_a*p_b*p_c*x*dsin(p_c*y)/t1
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_c**2*dlog(t1)*dcos(p_c*y)
            case(47)
               ex=dexp(p_b*x**2)
               ty=dtanh(p_c*y)
               t1=1.d0-ty**2
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*ty*ex*(2.d0*p_b+4.d0*p_b**2*x**2)
               if(ii.eq.2.and.jj.eq.1) &
                  res=2.d0*p_a*p_b*p_c*x*ex*t1
               if(ii.eq.2.and.jj.eq.2) &
                  res=-2.d0*p_a*p_c**2*ex*ty*t1
            case(48)
               if(ii.eq.1.and.jj.eq.1) &
                  res=-p_a*p_b**2*dsin(p_b*x)*y**3
               if(ii.eq.2.and.jj.eq.1) &
                  res=3.d0*p_a*p_b*dcos(p_b*x)*y**2
               if(ii.eq.2.and.jj.eq.2) &
                  res=6.d0*p_a*dsin(p_b*x)*y
            case(49)
               if(ii.eq.1.and.jj.eq.1) &
                  res=12.d0*p_a*x**2*dcos(p_b*y)
               if(ii.eq.2.and.jj.eq.1) &
                  res=-4.d0*p_a*p_b*x**3*dsin(p_b*y)
               if(ii.eq.2.and.jj.eq.2) &
                  res=-p_a*p_b**2*x**4*dcos(p_b*y)
            case(50)
               tx=dtanh(p_b*x**3)
               t1=1.d0-tx**2
               ey=dexp(p_c*y)
               if(ii.eq.1.and.jj.eq.1) &
                  res=p_a*ey*(6.d0*p_b*x*t1-&
                     9.d0*p_b**2*x**4*tx*t1*2.d0)
               if(ii.eq.2.and.jj.eq.1) &
                  res=3.d0*p_a*p_b*p_c*x**2*t1*ey
               if(ii.eq.2.and.jj.eq.2) &
                  res=p_a*p_c**2*tx*ey
            case default
               res=0.d0
         end select
         !***
         end
      end
