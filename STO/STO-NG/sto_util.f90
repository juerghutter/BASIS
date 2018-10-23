! **************************************************************************************************
MODULE sto_util
   USE kinds
   USE mathconstants

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: fit, test

CONTAINS

! **************************************************************************************************
   SUBROUTINE fit(n,l,ivar,xv,cv,err)
      INTEGER, INTENT(IN)                     :: n, l, ivar
      REAL(dp), DIMENSION(:), INTENT(IN)      :: xv
      REAL(dp), DIMENSION(:), INTENT(OUT)     :: cv
      REAL(dp), INTENT(OUT)                   :: err

      LOGICAL                                 :: docalc
      REAL(dp), DIMENSION(:), ALLOCATABLE     :: sg, ax
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: gg, ai
      REAL(dp)                                :: fs, fg, scal
      INTEGER                                 :: m, ng, i, j, info

      docalc = .TRUE.
      ng = l+1
      m = ivar
      cv(1:m) = 0.0_dp
      err = HUGE(1._dp)
      IF(xv(1) < 0._dp) docalc = .FALSE.
      DO i=2,m
         IF(xv(i) < 1.10_dp*xv(i-1)) docalc = .FALSE.
      END DO

      IF(docalc) THEN
         ALLOCATE(sg(m), ax(m))
         fs = SQRT ( 2._dp**(2*n+1) / fac(2*n) )
         fg = SQRT ( SQRT(2._dp/pi) * 2._dp**(2*ng+1) / dfac(2*ng-1) )
         ax(1:m) = fg * xv(1:m)**(0.5_dp*ng+0.25_dp)
         CALL sgints(n,l,xv,sg,m)
         sg(1:m) = sg(1:m) * ax(1:m) * fs
         ALLOCATE(gg(m,m), ai(m,m))
         gg = 0.0_dp
         DO i=1,m
            DO j=1,i
               gg(i,j) = 0.5_dp*gamma1(l+1)*(xv(i)+xv(j))**(-l-1.5_dp)*ax(i)*ax(j)
               gg(j,i) = gg(i,j)
            END DO
         END DO
         ai(1:m,1:m)=gg(1:m,1:m)
         CALL DPOTRF("U",m,ai,m,info)
         cv(1:m) = sg(1:m)
         CALL DPOTRS("U",m,1,ai,m,cv,m,info)
         scal = SUM(MATMUL(gg(1:m,1:m),cv(1:m))*cv(1:m))
         cv(1:m) = cv(1:m)/SQRT(scal)

         err = 2._dp - 2._dp*SUM(cv(1:m)*sg(1:m))

         DEALLOCATE(sg, ax, gg)
      END IF

   END SUBROUTINE fit

   SUBROUTINE test(n,l,ivar,xv,cv)
      INTEGER, INTENT(IN)                  :: n, l, ivar
      REAL(dp), DIMENSION(:), INTENT(IN)   :: xv, cv

   END SUBROUTINE test

   SUBROUTINE sgints(n,l,a,sval,ma)
      INTEGER, INTENT(IN)  :: n, l
      REAL(dp), DIMENSION(:), INTENT(IN) :: a
      REAL(dp), DIMENSION(:), INTENT(OUT) :: sval
      INTEGER, INTENT(IN)  :: ma

      INTEGER              :: i, m
      REAL(dp)             :: f, t, cost, sint, sint2, x, w
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: r,r2,wr,ff

      m = 1000
      f = pi/REAL(m+1, dp)
      ALLOCATE(r(m),r2(m),ff(m),wr(m))
      DO i = 1, m
         t = REAL(i, dp)*f
         cost = COS(t)
         sint = SIN(t)
         sint2 = sint**2
         x = REAL(2*i-m-1, dp)/REAL(m+1, dp)- &
                2.0_dp*(1.0_dp+2.0_dp*sint2/3.0_dp)*cost*sint/pi
         w = 16.0_dp*sint2**2/REAL(3*(m+1), dp)
         r(m+1-i) = LOG(2.0_dp/(1.0_dp-x))/LOG(2.0_dp)
         r2(m+1-i) = r(m+1-i)**2
         wr(m+1-i) = w*r2(m+1-i)/(LOG(2.0_dp)*(1.0_dp-x))
      END DO

      ff(1:m) = EXP(-r(1:m)) * r(1:m)**(n+l-1) * wr(1:m)
      DO i=1,ma
         sval(i) = SUM(ff(1:m)*EXP(-a(i)*r(1:m)**2))
      END DO

      DEALLOCATE(r,ff,wr)

   END SUBROUTINE sgints


END MODULE sto_util

