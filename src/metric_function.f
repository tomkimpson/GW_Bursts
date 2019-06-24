MODULE metric_function

  
  
  
USE SHARED_CONSTANTS
IMPLICIT NONE

!  Define access to subroutines.

PUBLIC :: metric_covar, metric_contra

CONTAINS

SUBROUTINE metric_covar(r, theta, &
                         g00,g01,g02,g03,&
                         g10,g11,g12,g13,&
                         g20,g21,g22,g23,&
                         g30,g31, g32,g33)

REAL(KIND=dp), INTENT(IN) :: r, theta
REAL(KIND=dp), INTENT(OUT) :: g00,g01,g02,g03,&
                              g10,g11,g12,g13,&
                              g20,g21,g22,g23,&
                              g30,g31,g32,g33
REAL(KIND=dp) :: Sg, Dl, M                         


M = 1.00_dp
Sg = r**2 + a**2*cos(theta)**2
Dl = r**2 + a**2 - 2.0_dp*r

g00 = -(1.0_dp - 2.0_dp*r/Sg)
g01 = 0.0_dp
g02 = 0.0_dp
g03 = -2.0_dp*r*a*sin(theta)**2/Sg

g11 = Sg/Dl
g12 = 0.0_dp
g13 = 0.0_dp

g22 = Sg
g23 = 0.0_dp

g33 = (r**2 + a**2 + 2.0_dp*r*a**2*sin(theta)**2/Sg)*sin(theta)**2

g10 = g01
g20 = g02
g30 = g03
g21 = g12
g31 = g13
g32 = g23

END SUBROUTINE metric_covar

SUBROUTINE metric_contra(r, theta, H00,H01,H02,H03,&
                                   H10,H11,H12,H13,&
                                   H20,H21,H22,H23,&
                                   H30,H31, H32,H33)
REAL(KIND=dp), INTENT(IN) :: r, theta
REAL(kind=dp), INTENT(OUT)::       H00,H01,H02,H03,&
                                   H10,H11,H12,H13,&
                                   H20,H21,H22,H23,&
                                   H30,H31, H32,H33
REAL(KIND=dp) :: Sg, Dl, M                          

M = 1.00_dp
Sg = r**2+a**2*cos(theta)**2  
Dl = r**2 + a**2 - 2.0_dp*M*r

H00 = -((r**2+a**2)+2.0_dp*M*r*a**2*sin(theta)**2/Sg)/Dl
H01 = 0.0_dp
H02 = 0.0_dp
H03 = -2.0_dp*M*r*a/(Sg*Dl)

H11 = Dl/Sg
H12 = 0.0_dp
H13 = 0.0_dp

H22 = 1.0_dp/Sg
H23 = 0.0_dp

H33 = (1.0_dp-2.0_dp*M*r/Sg)/(Dl*sin(theta)**2)


H10 = H01
H20 = H02
H30 = H03
H21 = H12
H31 = H13
H32 = H23

END SUBROUTINE metric_contra
END MODULE metric_function
