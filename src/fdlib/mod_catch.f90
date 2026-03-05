!***********************************************************************
      module mod_catch
!***********************************************************************
      implicit none
      public :: catch_error
      contains
      !*****************************************************************
         subroutine catch_error(err,msg,proc,comment)
      !*****************************************************************
         implicit none
         logical :: &
            err
         character(*) &
            msg,&
            proc
         character(*), optional :: &
            comment
         !***
         if(err) then
            write(*,'(/2x,a/)')    '!!! FATAL ERROR !!!'
            write(*,'(5x,a)')     'Message   : '//trim(msg)
            if(present(comment)) &
               write(*,'(5x,a)')  'Comment   : '//trim(comment)
            write(*,'(5x,a)')     'Procecure : '//trim(proc)
            write(*,*)
            stop
         endif
         !***
         end
      end
