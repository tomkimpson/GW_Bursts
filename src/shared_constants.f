MODULE SHARED_CONSTANTS
IMPLICIT NONE


!Define precision
integer, parameter :: dp = selected_real_kind(32,307)




!Universal constants
real(kind = dp), parameter :: PI = 4.D0*ATAN(1.D0)
real(kind=dp), PARAMETER :: light_c = 3d8
real(kind=dp), PARAMETER :: Newton_g = 6.67408d-11
real(kind=dp), PARAMETER :: Msun = 1.989d30
real(kind=dp), PARAMETER :: electron_charge = 4.80320425D-10 !CGS
real(kind=dp), PARAMETER :: electron_mass = 9.10938356D-28 !CGS
real(kind=dp), PARAMETER :: pc = 3.086e16 !Parsec in meters
real(kind=dp), PARAMETER :: au = 1.496d11 !AU in meters






!!!!---------System parameters which can be changed-------------!!!!!!!!!
!BH
real(kind=dp), PARAMETER :: MBH = 4.310D6!Solar masses
!real(kind=dp), PARAMETER :: MBH = 2.2310D3!Solar masses
real(kind=dp), PARAMETER :: a = 0.60_dp !BH spin parameter
real(kind=dp), PARAMETER :: Rhor = 1.0_dp+sqrt(1.0_dp-a**2) + 1d-2 !Horizon + eps
real(kind=dp), PARAMETER :: mu = Newton_g*MBH*Msun !m^3 s^-2
!PSR
real(kind=dp), PARAMETER :: MPSR = 1.4000_dp !Solar masses
real(kind=dp), PARAMETER :: RPSR = 1d4 !m
real(kind=dp), PARAMETER :: p0 = 1.00_dp !pulsar period in milliseconds


!Misc
real(kind=dp), PARAMETER :: convert_m = light_c**2/(Newton_g*MBH*Msun) 
real(kind=dp), PARAMETER :: convert_s = light_c**3/(Newton_g*MBH*Msun)
real(kind=dp), PARAMETER :: m0 = MPSR/MBH
real(kind=dp), parameter :: r0 = RPSR*convert_m !Pulsar radius in geo units


!Orbit
real(kind=dp), PARAMETER :: ecc = 0.90_dp





!IF SPECIFY RP

!real(kind=dp), parameter :: rp =15.60_dp
!real(kind=dp), parameter :: semi_latus = rp * (1.0_dp + ecc)
!real(kind=dp),PARAMETER :: ra = semi_latus/(1.0_dp - ecc)
!real(kind=dp) :: r_init = rp * 83.200_dp !For fig 1
!real(kind=dp), parameter :: semi_major = semi_latus/(1.0_dp - ecc**2.0_dp)


! IF SPECIFY ORBITAL PERIOD

real(kind=dp), PARAMETER :: periodKepler = 0.10_dp !Keplerian period in years
real(kind=dp), PARAMETER :: periodKeplerSec = periodKepler *365.0_dp*24.0_dp*3600.0_dp !in seconds
real(kind=dp), PARAMETER :: semi_major = convert_m * (periodKeplerSec**2 * Newton_g*MBH*Msun / (4.0_dp*PI**2.0_dp))**(1.0_dp/3.0_dp)


real(kind=dp),PARAMETER :: semi_latus = semi_major*(1.0_dp - ecc**2.0_dp)
real(kind=dp),PARAMETER :: ra = semi_latus/(1.0_dp - ecc)
real(kind=dp),PARAMETER :: rp = semi_latus/(1.0_dp + ecc)
real(kind=dp) :: r_init = rp !(ra+rp)/2.0_dp



real(kind=dp), PARAMETER :: BigTheta = 10.000_dp !This is i from Babak, Fang etc
real(kind=dp), PARAMETER :: theta_min = (90.0_dp - BigTheta) * PI/180.0_dp !Minimum latitude reached by PSR





real(kind=dp), PARAMETER :: sigma = 0.00_dp !orbital inclination w.r.t y-axis
real(kind=dp), PARAMETER  :: theta_init = sigma+PI/2.0_dp
real(kind=dp), PARAMETER :: phi_init = PI !0.0_dp !PI/2.0_dp !0.0_dp !4.0_dp*PI/4.0_dp !0.0_dp
!Components of rieman tensor





real(kind=dp), PARAMETER :: I = 0.40_dp*(MPSR*Msun)*RPSR**2 !moment of inertia
real(kind=dp), PARAMETER :: spin_factor = light_c/(Newton_g*(MBH*Msun)**2)
real(kind=dp), PARAMETER :: s0 = spin_factor*I*2*PI/(p0*1.0d-3) !S = I omega

INTEGER(kind=dp), PARAMETER :: circular = 0 !On/Off 1 = circular, 0= eccentric
real(kind=dp), PARAMETER :: xi = PI/4.0_dp, eta = PI/2.0_dp !0.0_dp !PI/2.0_dp 
!real(kind=dp), PARAMETER :: phi = PI/4.0_dp !0.0_dp !0.0_dp !Think of as (st,sx,sz,sy)
real(kind=dp), PARAMETER :: lambda = 1.00_DP
real(kind=dp), PARAMETER :: stheta = PI/ 4.0_dp 
real(kind=dp), PARAMETER :: sphi = PI/ 4.0_dp 



!Observational parameters
real(kind=dp), PARAMETER :: N_orbit = 1.50_dp !1.0_dp/periodKepler!10.0_dp !number of orbits




!!!!---------Integration constants-------------!!!!!!!!!

real(kind=dp), parameter :: integration_time = 1.0_dp !seconds
real(kind=dp), parameter :: fine_resolution =  1.0_dp !seconds
real(kind=dp), parameter :: coarse_resolution =  100.0_dp !seconds
real(kind=dp), parameter :: h_s = coarse_resolution
real(kind=dp), parameter :: escal = 1.0d20
real(kind=dp), parameter :: S = 0.90_dp
real(kind=dp), parameter :: Pgrow = -0.20_dp
real(kind=dp), parameter :: Pshrink = -0.250_dp
real(kind=dp), parameter :: errcon = (5.0_dp/S)**(1.0_dp/Pgrow)
real(kind=dp), parameter :: maxh = 0.50_dp * convert_s
real(kind=dp) :: h


!Cash-Karp parameters
real(kind = dp) :: B21=1.0/5.0
real(kind = dp) :: B31 = 3.0/40.0 , B32 = 9.0/40.0
real(kind = dp) :: B41 = 3.0/10.0, B42 = -9.0/10.0, B43 = 6.0/5.0 
real(kind = dp) :: B51 = -11.0/54.0, B52 = 5.0/2.0, B53 = -70.0/27.0
real(kind = dp) :: B54 = 35.0/27.0
real(kind = dp) :: B61 = 1631.0/55296.0, B62 = 175.0/512.0, B63 = 575.0/13824.0
real(kind = dp) :: B64 = 44275.0/110592.0, B65 = 253.0/4096.0
real(kind = dp) :: c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0
real(kind = dp) :: c6=512.0/1771.0
real(kind = dp) :: cbar1 = 2825.0/27648.0, cbar3 = 18575.0/48384.0
real(kind = dp) :: cbar4=13525.0/55296.0, cbar5 = 277.0/14336.0, cbar6 = 1.0/4.0


!Output 

character(len=200) :: PLOToutfile = '/unsafe/tok2/GravWaves/Trajectory.txt'
character(len=200) :: DATAoutfile = '/unsafe/tok2/GravWaves/Trajectory.dat'



!Observer Location - for setting the wave forms
real(kind=dp) :: OBSTheta != 0.0_dp !PI/2.0_dp !0.0_dp !PI/2.0_dp !p !PI/4.0_dp !PI/2.0_dp
real(kind=dp) :: OBSPhi = 0.0_dp
!real(kind=dp) :: OBSR = 4.0d3 * pc * convert_m !Distance to BH in rg
!real(kind=dp) :: OBSR = 763.0d3 * pc * convert_m !Distance to BH in rg
real(kind=dp) :: OBSR = 8.330d3 * pc * convert_m !Distance to BH in rg
real(kind=dp), dimension(3) :: Nvector
real(kind=dp),parameter :: Tobs = 11.570_dp * 24.0_dp * 3600.0_dp * convert_s
!Observation time. first numebr is days
!
!!!!!!!!!!!!!!!!!!-----------GLOBALS-----------!!!!!!!!!!!!!!!!!

real(kind=dp) :: FundFreq,timescale,vc

real(kind=dp) :: theta  !Think of as (st,sx,sz,sy)
!Globals
real(kind=dp) :: E, L, Q
REAL(KIND=dp) ::                                                                &
    R_0101, R_0102, R_0103, R_0112, R_0113, R_0123,                             &
    R_0201, R_0202, R_0203, R_0212, R_0213, R_0223,                             &
    R_0301, R_0302, R_0303, R_0312, R_0313, R_0323,                             &    
    R_1201, R_1202, R_1203, R_1212, R_1213, R_1223,                             &
    R_1301, R_1302, R_1303, R_1312, R_1313, R_1323,                             &
    R_2301, R_2302, R_2303, R_2312, R_2313, R_2323,                             &
! 
    R0_001, R0_002, R0_003, R0_012, R0_013, R0_023,                             &
    R0_101, R0_102, R0_103, R0_112, R0_113, R0_123,                             &
    R0_201, R0_202, R0_203, R0_212, R0_213, R0_223,                             &
    R0_301, R0_302, R0_303, R0_312, R0_313, R0_323,                             &
!
    R1_001, R1_002, R1_003, R1_012, R1_013, R1_023,                             &
    R1_101, R1_102, R1_103, R1_112, R1_113, R1_123,                             &
    R1_201, R1_202, R1_203, R1_212, R1_213, R1_223,                             &
    R1_301, R1_302, R1_303, R1_312, R1_313, R1_323,                             &
!
    R2_001, R2_002, R2_003, R2_012, R2_013, R2_023,                             &
    R2_101, R2_102, R2_103, R2_112, R2_113, R2_123,                             &
    R2_201, R2_202, R2_203, R2_212, R2_213, R2_223,                             &
    R2_301, R2_302, R2_303, R2_312, R2_313, R2_323,                             &
!
    R3_001, R3_002, R3_003, R3_012, R3_013, R3_023,                             &
    R3_101, R3_102, R3_103, R3_112, R3_113, R3_123,                             &
    R3_201, R3_202, R3_203, R3_212, R3_213, R3_223,                             &
    R3_301, R3_302, R3_303, R3_312, R3_313, R3_323



real(kind=dp) :: S_01, S_02, S_03, S_12, S_13, S_23

real(kind=dp) ::    G0_00, G0_01, G0_02, G0_03, G0_11, G0_12, G0_13, G0_22, G0_23, G0_33,       &
                    G1_00, G1_01, G1_02, G1_03, G1_11, G1_12, G1_13, G1_22, G1_23, G1_33,       &
                    G2_00, G2_01, G2_02, G2_03, G2_11, G2_12, G2_13, G2_22, G2_23, G2_33,       &
                    G3_00, G3_01, G3_02, G3_03, G3_11, G3_12, G3_13, G3_22, G3_23, G3_33
real(kind=dp) :: delta, RPS_0, RPS_1, RPS_2, RPS_3
real(kind=dp) :: start_factor
integer(kind=dp) :: factor, Printeger




END MODULE SHARED_CONSTANTS

