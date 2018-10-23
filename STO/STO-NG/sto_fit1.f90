Program sto
   USE kinds
   USE powell
   USE sto_util

   IMPLICIT NONE
 
   REAL(dp), DIMENSION(50) :: xv, cv, xx, xout
   REAL(dp)     :: xfac, x1, c1
   INTEGER      :: i,ivar,n,l,ic,ibeg,iend
   CHARACTER(len=1) :: cl
   TYPE(opt_state_type)                               :: ostate

   WRITE(6,FMT='(A)',ADVANCE="no") "STO n l "
   READ(5,*) n,cl 
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

   WRITE(6,FMT='(A)',ADVANCE="no") "Expansion ibeg iend " 
   READ(5,*) ibeg,iend

   DO ivar=ibeg,iend
      WRITE(6,*) " Expansion ",ivar

      ! initial guess
      IF(ibeg==iend) THEN
         WRITE(6,FMT='(A)',ADVANCE="no") "Geo Series x1 c1 " 
         READ(5,*) x1,c1
      ELSE
         x1 = 0.01_dp + (6-ivar)*0.01_dp
         c1 = 2.5_dp + 0.4_dp*(6-ivar)
      END IF
      xx(1) = x1
      xx(2) = c1
      ostate%nvar = 2
      ostate%rhoend = 1.e-10_dp
      ostate%maxfun = 500
      ostate%iprint = 1
      ostate%unit = 6
      ostate%nf = 0
      ostate%state = 0
      ostate%rhobeg = 0.01_dp 
      DO
         IF (ostate%state == 2) THEN
             CALL getgeo(xv,xx(1),xx(2),ivar)
             CALL fit(n,l,ivar,xv,cv,ostate%f)
         END IF
         IF (ostate%state == -1) EXIT
         CALL powell_optimize(ostate%nvar, xx, ostate)
      END DO
      ostate%state = 8
      CALL powell_optimize(ostate%nvar, xx, ostate)
      CALL getgeo(xv,xx(1),xx(2),ivar)
      CALL fit(n,l,ivar,xv,cv,ostate%f)
      WRITE(6,*) " Final Error Geometrical series ",ostate%f
      WRITE(6,*) " x=",xx(1),"     c=",xx(2)

      cv = 0.0_dp
      xv = 0.0_dp
      CALL getgeo(xv,xx(1),xx(2),ivar)
!!!!
!     xv(1) =   0.184830928571E-01_dp
!     xv(2) =   0.291824786372E-01_dp
!     xv(3) =   0.517397186484E-01_dp
!     xv(4) =   0.569412792726E-01_dp
!     xv(5) =   0.124987128687E+00_dp
!     xv(1) = 0.01_dp + (6-ivar)*0.01_dp
!     xfac = 2.5_dp + 0.5_dp*(6-ivar)
!     DO i=2,ivar
!        xv(i) = xfac * xv(i-1)
!     END DO
      CALL putexp(xv,xx,ivar)
      ostate%nvar = ivar
      ostate%rhoend = 1.e-14_dp
      ostate%maxfun = 5000
      ostate%iprint = 1
      ostate%unit = 6

      ostate%nf = 0
      ostate%state = 0
      ostate%rhobeg = 0.10_dp 

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
      CALL fit(n,l,ivar,xv,cv,ostate%f)
      WRITE(6,*) " Final Error ",ostate%f,"        Nfun:",ostate%nf

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
      xv(1) = xx(1)**2
      DO i=2,m
         xv(i) = xv(i-1) + xx(i)**2
      END DO
!     xv(1:m) = xx(1:m)
   END SUBROUTINE getexp

   SUBROUTINE putexp(xv,xx,m)
      REAL(dp),DIMENSION(:) :: xv, xx
      INTEGER               :: m, i
      xx(1) = SQRT(xv(1))
      DO i=2,m
         xx(i) = SQRT( xv(i) - xv(i-1) )
      END DO
!     xx(1:m) = xv(1:m)
   END SUBROUTINE putexp

   SUBROUTINE getgeo(xv,x1,c1,m)
      REAL(dp),DIMENSION(:) :: xv
      REAL(dp)              :: x1, c1
      INTEGER               :: m, i
      xv(1) = x1
      DO i=2,m
         xv(i) = xv(i-1) * c1
      END DO
   END SUBROUTINE getgeo

END Program sto
