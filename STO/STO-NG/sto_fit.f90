Program sto
   USE kinds
   USE powell
   USE sto_util

   IMPLICIT NONE
 
   REAL(dp), DIMENSION(50) :: xv, cv, xx, xout
   REAL(dp)     :: xfac, x1, c1, errin, errout
   INTEGER      :: i,ivar,n,l,ic,ibeg,iend
   CHARACTER(len=1) :: cl
   TYPE(opt_state_type)                               :: ostate

   READ(5,*) n,cl,ibeg,iend
   SELECT CASE (cl)
   CASE ("S","s")
      l = 0
   CASE ("P","p")
      l = 1
   CASE ("D","d")
      l = 2
   CASE ("F","f")
      l = 3
   CASE ("G","g")
      l = 4
   CASE ("H","h")
      l = 5
   CASE ("I","i","J","j")
      l = 6
   CASE ("K","k")
      l = 7
   CASE ("L","l")
      l = 8
   CASE ("M","m")
      l = 9
   CASE DEFAULT
      stop ("l val")
   END SELECT

   WRITE(6,'(i2,A1)') n, cl
   DO ivar=ibeg,iend
      WRITE(6,*) " Expansion ",ivar
      READ(5,*) 
      xv = 0.0_dp
      DO i=1,ivar
         READ(5,*) xv(i)
      END DO
      CALL order(xv,ivar)
      cv = 0.0_dp
      CALL fit(n,l,ivar,xv,cv,errin)
      WRITE(6,*) " Initial Error ",errin
      CALL putexp(xv,xx,ivar)
      ostate%nvar = ivar
      ostate%rhoend = 1.e-14_dp
      ostate%maxfun = 5000
      ostate%iprint = 1
      ostate%unit = 6

      ostate%f = 0.0_dp
      ostate%nf = 0
      ostate%state = 0
      ostate%rhobeg = 0.001_dp 

      DO

         IF (ostate%state == 2) THEN
             CALL getexp(xv,xx,ivar)
             CALL fit(n,l,ivar,xv,cv,ostate%f)
         END IF

         IF (ostate%state == -1) EXIT

         CALL powell_optimize(ostate%nvar, xx, ostate)

      END DO

      ostate%state = 8
      CALL powell_optimize(ostate%nvar, xx, ostate)
      CALL getexp(xv,xx,ivar)
      CALL fit(n,l,ivar,xv,cv,errout)
      WRITE(6,*) " Final Error   ",errout,"        Nfun:",ostate%nf

      DO i=1,ivar
         write(6,'(A,i1,A,E20.12,A,A,i1,A,E20.12,A)') "            alpha(",i,") = ",xv(i),"_dp",&
           "; coef(",i,") = ",cv(i),"_dp"
      END DO

      ! CALL test(n,l,ivar,xv,cv)

   END DO

   CONTAINS
   
   SUBROUTINE getexp(xv,xx,m)
      REAL(dp),DIMENSION(:) :: xv, xx
      INTEGER               :: m, i
!     xv(1) = xx(1)**2
!     DO i=2,m
!        xv(i) = xv(i-1) + xx(i)**2
!     END DO
      xv(1:m) = xx(1:m)
   END SUBROUTINE getexp

   SUBROUTINE putexp(xv,xx,m)
      REAL(dp),DIMENSION(:) :: xv, xx
      INTEGER               :: m, i
!     xx(1) = SQRT(xv(1))
!     DO i=2,m
!        xx(i) = SQRT( xv(i) - xv(i-1) )
!     END DO
      xx(1:m) = xv(1:m)
   END SUBROUTINE putexp

   SUBROUTINE order(xv,m)
      REAL(dp),DIMENSION(:) :: xv
      INTEGER               :: m, j
      REAL(dp),DIMENSION(:),ALLOCATABLE :: yy

      ALLOCATE(yy(m))
      yy = 0.0_dp
      DO i=1,m
         j=SUM(MINLOC(xv(1:m)))
         yy(i) = xv(j)
         xv(j) = HUGE(1._dp)
      END DO
      xv(1:m) = yy(1:m)
      DEALLOCATE(yy)

   END SUBROUTINE order

END Program sto
