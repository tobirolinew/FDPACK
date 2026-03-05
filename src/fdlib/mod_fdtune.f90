!***********************************************************************
      module mod_fdtune  !noise/offset estimator for finite differences
!***********************************************************************
      use iso_fortran_env, only : r64=>real64
      use mod_fdbase, only : &
        fdconfig, cdgrad, fwgrad, cdhess
      implicit none
      private
      public :: get_noise_table, get_fdoffsets
      contains
      !*****************************************************************
         function get_noise_table(&
                     func,v,maxoffsets,fixed,f,mnoise_scans,&
                     rtol_fround,rtol_offnoise,rtol_fnoise,&
                     rtol_fchange,cnvgd) &
                        result(noise_table)
      !*****************************************************************
      !
      !  This function estimates the noise level of a multivariate
      !  function for each free variable.
      !  It probes the function at increasing offsets and extrapolates
      !  the noise floor via Lagrange interpolation.
      !
      !  Returns:
      !    noise_table(:,1) - noise-floor offsets per variable
      !    noise_table(:,2) - estimated noise levels per variable
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
            !objective function:
         include 'src/iface/func_iface_inc.f90'
            !type for local variables:
         type lview_
            logical, allocatable :: fixed(:)
            logical :: cnvgd
            integer :: mnoise_scans=15
            real(r64) :: f
            real(r64) :: &
               rtol_fround=5.d0,&
               rtol_offnoise=0.1d0,&
               rtol_fnoise=10.d0,&
               rtol_fchange=10.d0
         end type
         type(lview_) :: &
            !lview_ instance:
               lview
         logical, optional :: &
            !flags for variable fixation:
               fixed(:),&
            !convergence flag:
               cnvgd
         integer, optional :: &
            !maximum number of noise scans:
               mnoise_scans
         real(r64), optional :: &
            !reference function value:
               f,&
            !relative tolerance for 
            !function rounding:
               rtol_fround,&
            !relative offset tolerance:
               rtol_offnoise,&
            !relative tolerance 
            !for function noise:
               rtol_fnoise,&
            !relative tolerance 
            !for function change:
               rtol_fchange
         real(r64) :: &
            !variable vector:
               v(:),&
            !maximum offsets:
               maxoffsets(:),&
            !roundoff-level 
            !function change:
               fround,&
            !current function change:
               fchange,&
            !trial function change:
               new_fchange,&
            !previous noise estimate:
               old_fnoise,&
            !previous function change:
               old_fchange,&
            !current probing offset:
               off,&
            !roundoff-level offset:
               roundoff,&
            !minimum safe probing offset:
               safeoff,&
            !trial probing offset:
               newoff,&
            !Lagrange extrapolation weight:
               weight
         real(r64), allocatable :: &
            !noise table containing 
            !the function noises
            !and the associated 
            !noise-floor offsets:
               noise_table(:,:),&
            !temporary variable vector:
               vtmp(:),&
            !sequence of sampled offsets:
               offsets(:),&
            !function changes 
            !at sampled offsets:
               fchanges(:)
         integer :: &
            !iterators:
               i,j,k,m,n
         !***
         associate(nv  => size(v),&
                   eps => epsilon(0.d0))
            !fill up the temporary 
            !variable vector
            vtmp=v

            !initialize the simple 
            !parameters in "lview"
            if(present(fixed)) then
               lview%fixed=fixed
            else
               lview%fixed=spread(.false.,1,nv)
            endif
            if(present(mnoise_scans)) &
               lview%mnoise_scans=mnoise_scans
            if(present(rtol_fround)) &
               lview%rtol_fround=rtol_fround
            if(present(rtol_offnoise)) &
               lview%rtol_offnoise=rtol_offnoise
            if(present(rtol_fnoise)) &
               lview%rtol_fnoise=rtol_fnoise
            if(present(rtol_fchange)) &
               lview%rtol_fchange=rtol_fchange

            !evaluate or store the reference 
            !function value
            if(present(f)) then
               lview%f=f
            else
               lview%f=func(vtmp)
            endif

            !check for size incompatibility errors
            call catch_error(&
                     err=nv.ne.size(lview%fixed).or.&
                         nv.ne.size(maxoffsets),&
                     msg='size incompatibility detected.',&
                     proc='get_noise_table')

            !initialize the noise table (function
            !noises and noise-floor offsets)
            allocate(noise_table(nv,2),source=0.d0)

            !loop over variables
            lview%cnvgd=.true.
            do i=1,nv
               !skip variable if fixed
               if(lview%fixed(i)) cycle

               !set up aliases for variable i
               associate(vi       => v(i),&
                         vti      => vtmp(i),&
                         noiseoff => noise_table(i,1),&
                         fnoise   => noise_table(i,2),&
                         maxoff   => maxoffsets(i))
                  !determine roundoff thresholds
                  !for function and variable changes
                  fround=eps*max(1.d0,dabs(lview%f))
                  roundoff=eps*max(1.d0,dabs(vi))

                  !minimum offset safely above roundoff
                  safeoff=roundoff/dsqrt(eps)

                  !start from the maximum probing offset
                  off=maxoff

                  !adjust offset to a measurable function change
                  !(balance roundoff and signal levels)
                  do n=1,lview%mnoise_scans
                     !probe function change at current offset
                     vti=vi+off
                     fchange=max(dabs(func(vtmp)-lview%f),fround)
                     vti=vi

                     !estimate new offset from signal magnitude
                     newoff=fround*off/fchange

                     !enforce lower bounds on offset
                     if(off.eq.maxoff) then
                        newoff=max(newoff,safeoff)
                     else
                        newoff=max(newoff,roundoff)
                     endif

                     !test new offset
                     vti=vi+newoff
                     new_fchange=&
                        max(dabs(func(vtmp)-lview%f),fround)
                     vti=vi

                     !accept offset if signal sufficient
                     !or relative change converged
                     if(new_fchange.ge.&
                           lview%rtol_fround*fround.or.&
                        dabs(newoff-off).lt.&
                           lview%rtol_offnoise*off) then
                        off=newoff
                        fchange=new_fchange
                        exit
                     endif

                     !geometric refinement of offset
                     off=dsqrt(off*newoff)
                  end do

                  !initialize local variables
                  !for the fine noise scan
                  offsets=[real(r64)::]
                  fchanges=[real(r64)::]

                  !initialize convergence monitors
                  old_fnoise=huge(0.d0)
                  old_fchange=huge(0.d0)

                  !scan increasing offsets and
                  !extrapolate the noise level
                  do m=n+1,lview%mnoise_scans
                     !store current scan level
                     offsets=[offsets,off]
                     fchanges=[fchanges,fchange]

                     !compute extrapolated noise estimate
                     fnoise=0.d0
                     do j=1,m-n
                        associate(offj => offsets(j),&
                                  fchj => fchanges(j))
                           weight=1.d0
                           do k=1,m-n
                              if(j.eq.k) cycle
                              associate(offk=>offsets(k))
                                 weight=weight*offk/(offj-offk)
                              end associate
                           end do
                           fnoise=fnoise+weight*fchj
                        end associate
                     end do

                     !bound noise estimate between
                     !roundoff level and measured change
                     fnoise=&
                        max(min(dabs(fnoise),fchange),fround)

                     !increase offset to next scan level
                     off=2.d0*lview%rtol_fnoise*off

                     !measure function change at new offset
                     vti=vi+off
                     fchange=&
                        max(dabs(func(vtmp)-lview%f),fround)
                     vti=vi

                     !check noise convergence
                     if(max(old_fnoise,fnoise).lt.&
                           lview%rtol_fnoise*&
                              min(fnoise,old_fnoise)) exit

                     !check function-change convergence
                     if(max(old_fchange,fchange).lt.&
                           lview%rtol_fchange*&
                              min(old_fchange,fchange)) exit

                     !update convergence monitors
                     old_fnoise=fnoise
                     old_fchange=fchange
                  end do

                  !store the minimum ("noise-floor") offset
                  if(m.eq.n+1) then
                     noiseoff=offsets(m-n)
                  else
                     noiseoff=offsets(m-n-1)
                  endif

                  !check convergence
                  lview%cnvgd=&
                     lview%cnvgd.and.&
                     m.le.lview%mnoise_scans
               end associate
            end do

            !save convergence flag if needed
            if(present(cnvgd)) cnvgd=lview%cnvgd
         end associate
         !***
         end

      !*****************************************************************
         function get_fdoffsets(&
                     func,v,fnoises,minoffsets,maxoffsets,fixed,f,&
                     approx,rtol_off,moff_scans,cnvgd) &
                        result(fdoffsets)
      !*****************************************************************
      !
      !  This function computes optimal finite-difference offsets for
      !  gradient and Hessian approximations of a noisy multivariate
      !  function.
      !  For each free variable (pair), it balances truncation error
      !  against noise amplification and iteratively refines a
      !  tight/loose bracket until the offset converges.
      !
      !  Returns: fdoffsets
      !    (i,i,1) - gradient offset for variable i
      !    (i,i,2) - diagonal Hessian offset for variable i
      !    (i,j,1) - i-direction offset for Hessian entry (i,j)
      !    (i,j,2) - j-direction offset for Hessian entry (i,j)
      !
      !*****************************************************************
         use mod_catch, only : catch_error
         implicit none
         include 'src/iface/func_iface_inc.f90'
            !type for local variables:
         type lview_
            logical, allocatable :: &
               fixed(:)
            real(r64) :: &
               f
            logical :: &
               cnvgd=.true.,&
               approx=.false.
            integer :: &
               moff_scans=1
            real(r64) :: &
               rtol_off=0.01d0
         end type
         type(lview_) :: &
            !lview_ instance:
               lview
         logical, optional :: &
            !fixed variable flags:
               fixed(:),&
            !approximate Hessian mode flag:
               approx,&
            !convergence flag:
               cnvgd
         integer, optional :: &
            !maximum number of offset scans:
               moff_scans
         real(r64), optional :: &
            !reference function value:
               f,&
            !relative offset convergence 
            !tolerance:
               rtol_off
         real(r64) :: &
            !variable vector:
               v(:),&
            !noise levels per variable:
               fnoises(:),&
            !lower bounds on probing offsets:
               minoffsets(:),&
            !upper bounds on probing offsets:
               maxoffsets(:)
         real(r64), allocatable :: &
            !optimal finite-difference 
            !offsets (output):
               fdoffsets(:,:,:),&
            !working copy of offsets used 
            !during evaluation:
               offsets(:,:,:),&
            !raw truncation estimates 
            !per variable (gradient pass):
               raw_truns(:),&
            !diagonal Hessian buffer 
            !(temporary):
               hdtmp(:),&
            !tight-bracket diagonal 
            !Hessian estimates:
               hdtight(:),&
            !loose-bracket diagonal 
            !Hessian estimates:
               hdloose(:)
         real(r64) :: &
            !tight- and loose-bracket offsets:
               offtight,offloose,&
            !derivative estimates 
            !at bracket endpoints:
               dertight,derloose,&
            !lower and upper offset bounds 
            !for current entry:
               minoff,maxoff
         integer :: &
            !finite-difference order 
            !(1=forward, 2=central):
               ord,&
            !number of offset scans 
            !for current entry:
               mscans
         !***
         call init()
         call find()
         !***
         contains
         !**************************************************************
            subroutine init()
         !**************************************************************
         !
         !  This subroutine initializes options, validates array sizes, 
         !  pre-allocates working arrays, and determines the 
         !  finite-difference order.
         !
         !**************************************************************
            implicit none
            !***
            associate(nv       => size(v),&
                      fdscheme => fdconfig%fdscheme)
               !initialize the simple 
               !parameters in "lview"
               if(present(fixed)) then
                  lview%fixed=fixed
               else
                  lview%fixed=spread(.false.,1,nv)
               endif
               if(present(approx)) &
                  lview%approx=approx
               if(present(rtol_off)) &
                  lview%rtol_off=rtol_off
               if(present(moff_scans)) &
                  lview%moff_scans=moff_scans

               !evaluate or store the reference 
               !function value
               if(present(f)) then
                  lview%f=f
               else
                  lview%f=func(v)
               endif

               !check for size incompatibility errors
               call catch_error(&
                       err=nv.ne.size(lview%fixed).or.&
                           nv.ne.size(minoffsets).or.&
                           nv.ne.size(maxoffsets).or.&
                           nv.ne.size(fnoises),&
                       msg='size incompatibility detected.',&
                       proc='get_fdoffsets::init')

               !pre-allocate working arrays
               allocate(offsets(nv,nv,2),source=0.d0)
               fdoffsets=offsets
               raw_truns=spread(0.d0,1,nv)
               hdtmp=spread(0.d0,1,nv)
               hdtight=spread(0.d0,1,nv)
               hdloose=spread(0.d0,1,nv)

               !determine finite-difference 
               !order from fdscheme
               select case(fdscheme)
                  case('forward')
                     ord=1
                  case('central')
                     ord=2
                  case default
                     call catch_error(&
                             err=.true.,&
                             msg='unknown fdscheme.',&
                             proc='get_fdoffsets::init')
               end select

               !approx mode is only available
               !for central differences
               call catch_error(&
                       err=lview%approx.and.ord.eq.1,&
                       msg='approx mode is not available '//&
                           'for forward differences.',&
                       proc='get_fdoffsets::init')
            end associate
            !***
            end

         !**************************************************************
            subroutine find()
         !**************************************************************
         !
         !  This subroutine finds gradient offsets for all free 
         !  variables, then Hessian offsets for all free variable pairs.
         !  For forward differences, Hessian offsets are inherited
         !  directly from the optimal gradient offsets.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable indices:
                  i,j
            !***
            associate(nv=>size(v))
               !find gradient offsets 
               !for all free variables
               do i=1,nv
                  if(lview%fixed(i)) cycle
                  call offscan(i)
               end do

               !set Hessian offsets depending 
               !on the scheme
               if(ord.eq.1) then
                  !inherit gradient offsets 
                  !for forward Hessian
                  do i=1,nv
                     if(lview%fixed(i)) cycle
                     do j=1,i
                        if(lview%fixed(j)) cycle
                        fdoffsets(i,j,:)=&
                           [fdoffsets(i,i,1),fdoffsets(j,j,1)]
                     end do
                  end do
               else
                  !optimize Hessian offsets 
                  !for all free variable pairs
                  do i=1,nv
                     if(lview%fixed(i)) cycle
                     do j=1,i
                        if(lview%fixed(j)) cycle
                        call offscan(i,j)
                     end do
                  end do
               endif

               !save convergence flag 
               !if requested
               if(present(cnvgd)) &
                  cnvgd=lview%cnvgd
            end associate
            !***
            end

         !**************************************************************
            subroutine offscan(i,j)
         !**************************************************************
         !
         !  This subroutine performs an iterative offset scan for a 
         !  single gradient (j absent) or Hessian (j present) entry.
         !  It computes the optimal offset and refines the bracket
         !  until convergence or the scan limit is reached.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i,&
               !scan counter:
                  iscan
            integer, optional :: &
               !second variable index 
               !(Hessian only):
                  j
            real(r64) :: &
               !current candidate offset:
                  off,&
               !derivative estimate 
               !at offset "off":
                  der
            !***
            call setup(i,j)
            iscan=1
            do
               off=eval_offset(i,j)

               !exit if offset has converged
               if(min(off,offtight).gt.lview%rtol_off*&
                  max(off,offtight)) exit

               !increment scan counter 
               !and check limit
               iscan=iscan+1
               if(iscan.eq.mscans+1) exit

               !refine bracket with updated 
               !derivative estimate
               call update_offset(i,j,off)
               der=eval_der(i,j)
               call update_offset(i,j)
               call update_bracket(off,der)
            end do

            !store the converged offset
            call store(off,i,j)

            !check convergence
            lview%cnvgd=&
               lview%cnvgd.and.&
                  iscan.le.lview%moff_scans
            !***
            end

         !**************************************************************
            subroutine setup(i,j)
         !**************************************************************
         !
         !  This subroutine initializes the tight/loose offset bracket
         !  and derivative estimates for a gradient (j absent) or
         !  Hessian (j present) entry.
         !  Diagonal Hessian entries reuse stored gradient bracket
         !  data; off-diagonal entries are evaluated or approximated.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index 
               !(Hessian only):
                  j
            !***
            if(present(j)) then
               !restrict the maximum number of scans
               mscans=merge(1,lview%moff_scans,lview%approx)

               !initialize bracket from geometric means
               associate(minoffj => minoffsets(j),&
                         maxoffj => maxoffsets(j))
                  offtight=dsqrt(minoffj*maxoffj)
                  offloose=dsqrt(offtight*maxoffj)
                  minoff=minoffj
                  maxoff=maxoffj
               end associate

               !initialize derivative estimates 
               !at bracket endpoints
               if(i.eq.j) then
                  !reuse diagonal gradient 
                  !bracket data
                  dertight=hdtight(i)
                  derloose=hdloose(i)
               elseif(lview%approx) then
                  !skip evaluation in 
                  !approximate mode
                  dertight=huge(0.d0)
                  derloose=huge(0.d0)
               else
                  !evaluate derivative
                  !at both bracket points
                  call update_offset(i,j,offtight)
                  dertight=eval_der(i,j)
                  call update_offset(i,j,offloose)
                  derloose=eval_der(i,j)
                  call update_offset(i,j)
               endif
            else
               associate(minoffi => minoffsets(i),&
                         maxoffi => maxoffsets(i))
                  mscans=lview%moff_scans
                  !initialize bracket 
                  !from geometric means
                  offtight=dsqrt(minoffi*maxoffi)
                  offloose=dsqrt(offtight*maxoffi)
                  minoff=minoffi
                  maxoff=maxoffi
               end associate

               !evaluate derivative and store Hessian 
               !diagonal at both bracket endpoints
               call update_offset(i,off=offtight)
               dertight=eval_der(i)
               hdtight(i)=hdtmp(i)
               call update_offset(i,off=offloose)
               derloose=eval_der(i)
               hdloose(i)=hdtmp(i)
               call update_offset(i)
            endif
            !***
            end

         !**************************************************************
            subroutine store(off,i,j)
         !**************************************************************
         !
         !  This subroutine saves the converged offset into "fdoffsets".
         !  For gradient entries (j absent), it also stores the raw
         !  truncation estimate for use in approximate Hessian mode.
         !  For Hessian entries (j present), it saves both the i- and
         !  the j-direction offsets for the gradient.
         !
         !**************************************************************
            implicit none
            real(r64) :: &
               !converged offset:
                  off
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable 
               !index (Hessian only):
                  j
            !***
            associate(offi=>fdoffsets(i,i,1))
               if(present(j)) then
                  fdoffsets(i,j,:)=[offi,off]
               else
                  offi=off
                  raw_truns(i)=eval_raw_trun(i)
               endif
            end associate
            !***
            end

         !**************************************************************
            function eval_offset(i,j) result(omega)
         !**************************************************************
         !
         !  This function computes the optimal offset by clamping the
         !  raw truncation coefficient via the switch function and
         !  applying the unified formula
         !  omega = (S_pq(kappa,nu,1)/trun)^(1/(p+q)).
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index
               !(Hessian only):
                  j
            real(r64) :: &
               !optimal offset:
                  omega,&
               !raw truncation coefficient:
                  raw_trun,&
               !lower bound on truncation
               !coefficient:
                  min_trun,&
               !upper bound on truncation
               !coefficient:
                  max_trun,&
               !clamped truncation
               !coefficient:
                  trun,&
               !switch function
               !numerator (off=1):
                  numer
            integer :: &
               !noise exponent:
                  nexp
            !***
            !determine the noise exponent
            if(present(j)) then
               nexp=merge(1,2,i.ne.j)
            else
               nexp=1
            endif

            !evaluate the truncation
            !coefficient
            raw_trun=eval_raw_trun(i,j)

            !evaluate switch function
            !at key offsets
            numer=eval_switch(1.d0,i,j)
            min_trun=eval_switch(maxoff,i,j)
            max_trun=eval_switch(minoff,i,j)

            !clamp truncation coefficient
            !via switch function
            trun=max(min_trun,min(raw_trun,max_trun))

            !compute optimal offset
            omega=(numer/trun)**(1.d0/(ord+nexp))
            !***
            end            

         !**************************************************************
            function eval_der(i,j) result(der)
         !**************************************************************
         !
         !  This function returns a gradient component (j absent) or
         !  Hessian entry (j present) using the central difference
         !  scheme (forward differences do not optimize Hessian 
         !  offsets).
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index 
               !(Hessian only):
                  j
            real(r64), allocatable :: &
               !temporary gradient vector:
                  gtmp(:),&
               !temporary Hessian matrix:
                  htmp(:,:)
            real(r64) :: &
               !derivative estimate:
                  der
            !***
            if(present(j)) then
               !compute central Hessian and 
               !extract entry (i,j)
               if(ord.eq.2) then
                  htmp=cdhess(func,v,offsets,lview%f)
                  der=htmp(i,j)
               else
                  call catch_error(&
                          err=.true.,&
                          msg='Hessian offset optimization is '//&
                              'not available for forward '//&
                              'differences.',&
                          proc='get_fdoffsets::eval_der')
               endif
            else
               !compute gradient and extract
               !the ith gradien component
               if(ord.eq.2) then
                  gtmp=cdgrad(func,v,offsets,lview%f,hdtmp)
               else
                  gtmp=fwgrad(func,v,offsets,lview%f)
               endif
               der=gtmp(i)
            endif
            !***
            end

         !**************************************************************
            function eval_switch(off,i,j) result(switch)
         !**************************************************************
         !
         !  This function evaluates the switch function
         !  S_pq(kappa,nu,off) = q*kappa*nu / (p*off^(p+q)),
         !  retrieving the truncation coefficient at a given offset.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index
               !(Hessian only):
                  j
            real(r64) :: &
               !offset at which to evaluate
               !the switch function:
                  off,&
               !switch function value:
                  switch,&
               !stencil noise factor:
                  kappa,&
               !noise level:
                  nu
            integer :: &
               !noise exponent:
                  nexp
            !***
            associate(fni  => fnoises(i),&
                      offi => fdoffsets(i,i,1))
               !determine noise exponent
               if(present(j)) then
                  nexp=merge(2,1,i.eq.j)
               else
                  nexp=1
               endif

               !determine "kappa" and "nu"
               if(present(j)) then
                  associate(fnj=>fnoises(j))
                     if(i.eq.j) then
                        !pure Hessian (central only)
                        kappa=4.d0*dsqrt(6.d0)
                        nu=fni
                     else
                        !mixed Hessian (central only)
                        kappa=8.d0
                        nu=dsqrt(fni**2+fnj**2)/offi
                     endif
                  end associate
               else
                  !gradient
                  kappa=dsqrt(2.d0)
                  nu=fni
               endif

               !calculate the switch value
               switch=nexp*kappa*nu/(ord*off**(ord+nexp))
            end associate
            !***
            end            

         !**************************************************************
            function eval_raw_trun(i,j) result(rtrun)
         !**************************************************************
         !
         !  This function estimates the raw truncation coefficient for 
         !  a gradient (j absent) or Hessian (j present) entry.
         !  Richardson extrapolation is used from the bracket endpoints,
         !  or a gradient-based guess in approximate mode.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index 
               !(Hessian only):
                  j
            real(r64) :: &
               !raw truncation 
               !coefficient:
                  rtrun
            !***
            if(lview%approx.and.present(j)) then
               !gradient-based approximation
               !for off-diagonal entries
               associate(offi   => fdoffsets(i,i,1),&
                         offj   => fdoffsets(j,j,1),&
                         rtruni => raw_truns(i),&
                         rtrunj => raw_truns(j))
                  rtrun=min(rtruni/offj,rtrunj/offi)
               end associate
            else
               !Richardson extrapolation 
               !from bracket endpoints
               rtrun=dabs(derloose-dertight)/&
                          (offloose**ord-offtight**ord)
            endif
            !***
            end

         !**************************************************************
            subroutine update_offset(i,j,off)
         !**************************************************************
         !
         !  This subroutine updates the working offset array for a
         !  gradient (j absent) or Hessian (j present) entry.
         !  It resets to zero if "off" is absent.
         !
         !**************************************************************
            implicit none
            integer :: &
               !variable index:
                  i
            integer, optional :: &
               !second variable index 
               !(Hessian only):
                  j
            real(r64), optional :: &
               !new offset value 
               !(reset to zero if absent):
                  off
            !***
            if(present(j)) then
               if(present(off)) then
                  offsets(i,j,:)=[fdoffsets(i,i,1),off]
               else
                  offsets(i,j,:)=0.d0
               endif
            else
               if(present(off)) then
                  offsets(i,i,1)=off
               else
                  offsets(i,i,1)=0.d0
               endif
            endif
            !***
            end

         !**************************************************************
            subroutine update_bracket(off,der)
         !**************************************************************
         !
         !  This subroutine assigns the new offset and derivative to
         !  the nearer bracket endpoint, then re-sorts so that
         !  offtight <= offloose.
         !
         !**************************************************************
            use mod_swap, only : swap
            implicit none
            real(r64) :: &
               !new offset to insert 
               !into bracket:
                  off,&
               !derivative estimate 
               !at offset "off":
                  der
            !***
            !replace the nearer bracket endpoint
            if(dabs(off-offtight).gt.dabs(off-offloose)) then
               offtight=off
               dertight=der
            else
               offloose=off
               derloose=der
            endif

            !ensure offtight <= offloose
            if(offtight.gt.offloose) then
               call swap(offtight,offloose)
               call swap(dertight,derloose)
            endif
            !***
            end
         !***
         end
      !***
      end
