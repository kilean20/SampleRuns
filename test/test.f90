   Program DebugIt
   implicit none
   double precision:: q,identity
   q = 0.5d0
! Works:
   write(*,*) 1.0d0
   write(*,*) identity(q)
! Fails:
   write(*,*) 1.0d0,identity(q)
   End


   function identity(dum)
   implicit none
   double precision:: dum
   double precision:: identity
   write(*,*) 'Inside function.'
   identity = dum
   end function identity
