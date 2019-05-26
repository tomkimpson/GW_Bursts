PROGRAM kerr

USE SHARED_CONSTANTS
USE init_EQL
USE initial_conditions
USE runge_kutta 
IMPLICIT NONE
real(kind=dp) :: S0_init, S1_init, S2_init, S3_init
real(kind=dp) :: P0_init, P1_init, P2_init, P3_init
real(kind=dp) :: t_init, t_orbit,t_final
!Calculate initial conditions for circular motion. See Raine/Thomas
INTEGER, PARAMETER :: entries = 12, extra = 23
INTEGER(kind=dp) :: N_step
REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: yy
REAL(KIND=dp), DIMENSION(entries) :: Y_init
real(kind=dp) :: OmK, period,SMA_meters,SMA2
real(kind=dp) :: inclination
character(len=20) :: num_obs_str,stheta_str, theta_str
character(len= 100) ::eccentricity_str,p_str,I_str 
real(kind=dp) :: FinalPhi

print *, 'Periapsis approach:', rp
print *, 'Num Orbits:', N_orbit

!inclinatio print *, atan(sqrt(Q)/L) * 180.0_dp/PI









!Set initial stepsize
h = h_s*convert_s


!Calculate fundamental frequency
FundFreq = sqrt(((1.0_dp+ecc)*mu)/((rp/convert_m)**3.0_dp))

!Calculate fundamental timescale

vc = sqrt((1.0_dp + ecc)/rp)
timescale = (rp/vc)/convert_s !critical timescale in seconds


!Calculate initial EQL
if (circular .EQ. 1) then
  r_init = 8.0_dp
  call EQL_circular()
elseif (circular .EQ. 0) then  
  call EQL_eccentric(rp,ra,theta_min)
endif





!Calculate spatial components of spin in coordinate basis
if (stheta .EQ. 0.0_dp) then
S1_init = 0.0_dp
S2_init = 0.0_dp
S3_init = 0.0_dp
else
S1_init = s0*sin(stheta)*cos(sphi)
S2_init = -s0*cos(stheta)/r_init
S3_init = s0*sin(stheta)*sin(sphi)/(r_init*sin(theta_init))
endif





!Set initial conditions
print *, 'start IC'
call init_cond(S0_init,S1_init, S2_init, S3_init,P0_init, P1_init, P2_init,P3_init)

!Set up array
Y_init(1) = t_init
Y_init(2) = r_init
Y_init(3) = theta_init
Y_init(4) = phi_init
Y_init(5) = P0_init
Y_init(6) = P1_init
Y_init(7) = P2_init
Y_init(8) = P3_init
Y_init(9) = S0_init
Y_init(10) = S1_init
Y_init(11) = S2_init
Y_init(12) = S3_init

print *, 'START RK'

!Run the RK integrator
call rk(Y_init,yy,N_step,t_final)

print *, '------ MPD EXITED OK------'

END PROGRAM kerr
