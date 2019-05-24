MODULE check_module

USE SHARED_CONSTANTS
IMPLICIT NONE

CONTAINS

SUBROUTINE checks(RR,TT)
IMPLICIT NONE
real(kind = dp) :: RR, TT


if (RR .LT. 0 .and. abs(RR) < 1d-16) then
print *, ' RR is negative but small. Correction applied :', RR, 0.00_dp
RR = 0.00_dp
else if (RR .LT. 0 .and. abs(RR) > 1d-16) then
print *, 'RR is negative and big. Something not right :', RR
STOP
endif


if (TT .LT. 0 .and. abs(TT) < 1d-16) then
print *, ' TT is negative but small. Correction applied :', TT, 0.00_dp
TT = 0.00_dp
else if (TT .LT. 0 .and. abs(TT) > 1d-16) then
print *, 'TT is negative and big. Something not right :', TT
STOP
endif


END SUBROUTINE checks



SUBROUTINE check_dixon(r, S0, S1,S2,S3, &
                          P0,P1,P2,P3, &
                          g00,g01,g02,g03,&
                          g10,g11,g12,g13,&
                          g20,g21,g22,g23,&
                          g30,g31, g32,g33, &
                          m0)


IMPLICIT NONE
real(kind=dp) :: r, S0, S1, S2, S3, P0,P1,P2,P3
real(kind=dp) :: g00,g01,g02,g03,&
                 g10,g11,g12,g13,&
                 g20,g21,g22,g23,&
                 g30,g31, g32,g33
real(kind=dp) :: m_sq, s_sq,m0

if (r .LE. 1) then
    print *, 'The pulsar has crashed into the BH - try a smaller step size or else alter orbital parameters'
    STOP
endif


m_sq = -(g00*P0*P0 + g11*P1*P1 + g22*P2*P2 + g33*P3*P3              &
+ g01*(P0*P1 + P1*P0) + g02*(P0*P2 + P2*P0) + g03*(P0*P3 + P3*P0)   & 
+ g12*(P1*P2 + P2*P1) + g13*(P1*P3 + P3*P1) + g23*(P2*P3 + P3*P2))

s_sq = g00*S0*S0 + g11*S1*S1 + g22*S2*S2 + g33*S3*S3                &
+ g01*(S0*S1 + S1*S0) + g02*(S0*S2 + S2*S0) + g03*(S0*S3 + S3*S0)   &
+ g12*(S1*S2 + S2*S1) + g13*(S1*S3 + S3*S1) + g23*(S2*S3 + S3*S2)


if (m_sq .LT. 0) then
print *,m_sq, 'Error! The square of the mass is negative. See dixon_eqs.f and checks. Best to lower the step size'
STOP
endif
if (s_sq .LT. 0) then
print *, 'Error! The square of the spin is negative. See dixon_eqs.f and checks'
STOP
endif






m0 = sqrt(m_sq)

END SUBROUTINE check_dixon









END MODULE check_module
