MODULE initial_conditions


USE SHARED_CONSTANTS
USE check_module
USE metric_function
IMPLICIT NONE

CONTAINS

SUBROUTINE init_cond(S0, S1,S2,S3,Pt,P1,P2,P3)
IMPLICIT NONE
real(kind=dp) :: r,theta, sigma,delta, PP, RR,TT
real(kind=dp) :: S0, S1, S2, S3
real(kind=dp) :: tdot, rdot, thetadot, phidot
real(kind=dp) :: Pt, P1, P2,P3
real(kind=dp) :: m_sq, s_sq
REAL(KIND=dp) :: g00,g01,g02,g03,&
                 g10,g11,g12,g13,&
                 g20,g21,g22,g23,&
                 g30,g31,g32,g33
r = r_init
theta = theta_init

sigma = r**2+a**2*cos(theta)
delta= r**2+a**2-2.0_dp*r

PP = E*(r**2+a**2) - a*L
RR = PP**2-delta*(r**2 +Q+(L-a*E)**2)
TT = Q - cos(theta)**2 * (a**2*(1.00_dp-E**2)+L**2/sin(theta)**2)

call checks(RR,TT)

!Kerr diff equations - think of as magnitudes
tdot = (a*(L - a*E*sin(theta)**2) + (r**2 + a**2)*PP/delta)/sigma
rdot = sqrt(RR)/sigma !I HAVE CHANGED THIS TO A MINUS SIGN
thetadot = sqrt(TT)/sigma
phidot = ((L/sin(theta)**2 -a*E) + a*PP/delta)/sigma



!Momentum


if (circular .EQ. 1) then
Pt = m0*tdot
P1 = m0*rdot*sin(xi)*cos(eta)
P2 = -m0*thetadot*cos(xi)
P3 = m0*phidot*sin(xi)*sin(eta)
else if (circular .EQ. 0) then

Pt = m0*tdot
P1 = m0*rdot*cos(eta)
P2 = -m0*thetadot
P3 = m0*phidot*sin(eta)



endif






call  metric_covar(r, theta, &
                   g00,g01,g02,g03,&
                   g10,g11,g12,g13,&
     g20,g21,g22,g23,&
                   g30,g31, g32,g33)

!Zeroth Component of spin
S0 = -((g01*P0 + g11*P1 + g12*P2 + g13*P3)*S1  &
+ (g02*P0 + g12*P1 + g22*P2 + g23*P3)*S2       &
+ (g03*P0 + g13*P1 + g23*P2 + g33*P3)*S3)/(g00*P0 + g01*P1 + g02*P2 + g03*P3)








!Keplerian orbital frequency

!Calculate m^2 and s^2
m_sq = -(g00*P0**2 + g11*P1**2 + g22*P2**2 + g33*P3**2                                       &
+ 2.0_DP*(g01*P0*P1 + g02*P0*P2 + g03*P0*P3 + g12*P1*P2 + g13*P1*P3 + g23*P2*P3))

s_sq = (g00*S0**2 + g11*S1**2 + g22*S2**2 + g33*S3**2                                        &
+ 2.0_DP*(g01*S0*S1 + g02*S0*S2 + g03*S0*S3 + g12*S1*S2 + g13*S1*S3 + g23*S2*S3))

END SUBROUTINE init_cond


subroutine printstatus(period,hs,h)

real(kind=dp) :: period, hs,h
character(len=10) :: step_string, orbit_string

write(orbit_string, '(f10.2)') N_orbit
write(step_string, '(f10.2)') h_s


Print *, '------ START MPD RUN-----'
Print *, '------ Parameters--------'
print *, 'Semi major axis: ', semi_major, ' rg'
print *, 'Eccentricity: ', ecc
print *, 'Period', period/(convert_s*60*60*24*7), ' weeks'
print *, 'Initial stepsize = ', step_string, ' seconds'
print *, 'N orbits:', orbit_string

end subroutine printstatus






END MODULE initial_conditions
