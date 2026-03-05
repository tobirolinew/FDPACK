      interface
         subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         use iso_fortran_env, only : r64=>real64
         implicit none
         character :: jobz,uplo
         integer :: n, lda,lwork
         real(r64) :: a(lda,*),work(*)
         real(r64) :: w(*)
         integer :: info
         end
         subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,&
                           work,lwork,info)
         implicit none
         character :: jobu,jobvt
         integer :: m,n,lda,ldu,ldvt,lwork,info
         real*8 :: a(lda,*),s(*),u(ldu,*),vt(ldvt,*),work(*)
         end
         subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
         implicit none
         integer :: m,n,incx,incy,lda
         real*8 :: alpha,x(*),y(*),a(lda,*)
         end
      end interface
