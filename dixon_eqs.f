MODULE Dixon_equations

USE SHARED_CONSTANTS
USE tensors
USE metric_function
USE check_module
IMPLICIT NONE
PRIVATE

!  Define access to subroutines.

PUBLIC :: derivs

CONTAINS


SUBROUTINE derivs(y,dydx)
        
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y
REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: dydx


REAL(KIND=dp) :: t, r, theta, phi, &
                 P0, P1, P2, P3,&
                 S0, S1, S2, S3,&
                 m0, M                                        

real(kind=dp) :: g00,g01,g02,g03,&
                 g10,g11,g12,g13,&
                 g20,g21,g22,g23,&
                 g30,g31, g32,g33


real(kind=dp) :: H00,H01,H02,H03,&
                 H10,H11,H12,H13,&
                 H20,H21,H22,H23,&
                 H30,H31, H32,H33
                     
real(kind=dp) :: Pdot_0, Pdot_1, Pdot_2, Pdot_3

real(kind=dp) :: Sdot_0, Sdot_1, Sdot_2, Sdot_3


real(kind=dp) :: V0, V1, V2, V3

! Read in t, r, theta, phi, P0, P1, P2, P3, S0, S1, S2, S3
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

M = 1.00_dp


!Calculate covariant metric components (i.e. lower indices)
call metric_covar(r, theta, g00,g01,g02,g03,&
                            g10,g11,g12,g13,&
                            g20,g21,g22,g23,&
                            g30,g31, g32,g33)


call metric_contra(r, theta, H00,H01,H02,H03,&
                             H10,H11,H12,H13,&
                             H20,H21,H22,H23,&
                             H30,H31, H32,H33)
 

!Check everything is OK
call check_dixon(r, S0, S1,S2,S3, &
                    P0,P1,P2,P3, &
                    g00,g01,g02,g03,&
                    g10,g11,g12,g13,&
                    g20,g21,g22,g23,&
                    g30,g31, g32,g33, &
                    m0)


! Calculate components of antisymmetric spin tensor


call calc_spin_antisym(r,theta, &
                       S0,S1,S2,S3,&
                       P0,P1,P2,P3,&
                       H00,H01,H02,H03,&
                       H10,H11,H12,H13,&
                       H20,H21,H22,H23,& 
                       H30,H31, H32,H33)


!Calculate components of Christoffels symbols
call calc_christoffel(r,theta)

!Calculate components of Riemann tensor
call calc_riemann(r,theta)

!Calculate 4-velocity

call calc_FourVel(m0,P0,P1,P2,P3, &
                     V0,V1,V2,V3, &
                     g00,g11,g22,g33,g30)



!Calculate 4-momentum
call calc_FourMom(V0,V1,V2,V3, &
                  P0,P1,P2,P3, &
                  Pdot_0, Pdot_1,Pdot_2,Pdot_3)

!Calculate 4-spin


call calc_FourSpin(V0,V1,V2,V3,&
                         S0,S1,S2,S3,&
                         P0,P1,P2,P3,&
                         m0,&
                         Sdot_0, Sdot_1, Sdot_2, Sdot_3)






!print *, '4 velocity magnitude:', g00*V0*V0 + &
 !                                 g11*V1*V1 + &
  !                                g22*V2*V2 + &
   !                               g33*V3*V3 + &
    !                              2.0_dp*g30*V0*V3



!print *, '4 momentum magnitude:', g00*P0*P0 + &
 !                                 g11*P1*P1 + &
  !                                g22*P2*P2 + &
   !                               g33*P3*P3 + &
    !                              2.0_dp*g30*P0*P3



!Save to array
dydx(1) = V0
dydx(2) = V1
dydx(3) = V2
dydx(4) = V3
dydx(5) = Pdot_0
dydx(6) = Pdot_1
dydx(7) = Pdot_2
dydx(8) = Pdot_3
dydx(9) = Sdot_0
dydx(10) = Sdot_1
dydx(11) = Sdot_2
dydx(12) = Sdot_3

call metric_covar(r, theta, g00,g01,g02,g03,&
                            g10,g11,g12,g13,&
                            g20,g21,g22,g23,&
                            g30,g31, g32,g33)

END SUBROUTINE derivs






END MODULE Dixon_equations
