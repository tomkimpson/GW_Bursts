MODULE spin_condition
USE SHARED_CONSTANTS
USE metric_function
IMPLICIT NONE
PRIVATE

!  Define access to subroutines.

PUBLIC :: spin_cond

CONTAINS

SUBROUTINE spin_cond(a, y)

REAL(KIND=dp), INTENT(IN) :: a
REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: y
REAL(KIND=dp) :: t, r, theta, phi, P0, P1, P2, P3, S0, S1, S2, S3
REAL(KIND=dp) :: g00,g01,g02,g03,&
                 g10,g11,g12,g13,&
                 g20,g21,g22,g23,&
                 g30,g31,g32,g33



t = y(1)
r = y(2)
theta = y(3)
phi = y(4)
P0 = y(5)
P1 = y(6)
P2 = y(7)
P3 = y(8)
S0 = y(9)
S1 = y(10)
S2 = y(11)
S3 = y(12)


call  metric_covar(r, theta, &
                   g00,g01,g02,g03,&
                   g10,g11,g12,g13,&
                   g20,g21,g22,g23,&
                   g30,g31, g32,g33)

!  Overwrite S0 with value derived from spin condition.

S0 = -((g01*P0 + g11*P1 + g12*P2 + g13*P3)*S1  &
+ (g02*P0 + g12*P1 + g22*P2 + g23*P3)*S2       &
+ (g03*P0 + g13*P1 + g23*P2 + g33*P3)*S3)/(g00*P0 + g01*P1 + g02*P2 + g03*P3)


y(9) = S0

END SUBROUTINE spin_cond

END MODULE spin_condition
