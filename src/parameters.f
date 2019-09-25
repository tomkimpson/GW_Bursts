module parameters
implicit none


!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)
!integer, parameter :: dp = selected_real_kind(15,307)

!Define path to directory where you want the output files to go
!character(len=200) :: PathOut = '/Users/tomkimpson/Data/SpinCurv/'


!System Parameters - these can be changed to model different sorts of orbits

!The semi-major axis of the orbit can be defined either directly, or related to the period via Keplers 3rd

real(kind=dp), parameter :: KeplerianPeriod = 0.10_dp
real(kind=dp), parameter :: r_set = 6.0_dp !Set the sma. Only used if KepPer eq 0

real(kind=dp), parameter :: circular = 0 !Turn on/off a precessing circular orbit. If on subsequent parameters are ignored

real(kind=dp), parameter :: eccentricity = 0.90_dp !Orbital eccentricity
real(kind=dp), parameter :: iota = 15.0_dp !Inclination w.r.t equatorial plane in degrees
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar masses
real(kind=dp), parameter :: a=0.60_dp !BH spin parameter
real(kind=dp), parameter :: MPSR = 1.400_dp !pulsar mass in solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !pulsar radius in km
real(kind=dp), parameter :: p0 = 1.0e-3_dp !pulsar spin period in seconds
real(kind=dp), parameter :: N_orbit = 1.0_dp !Number of orbits to integrate
real(kind=dp), parameter :: lambda = 1.0_dp !Turn on/off spin-curvature coupling (1 = on)

!For GW calcs
real(kind=dp), parameter :: observer_distance = 8.330d3 !parsecs


!Format parameters
integer(kind=dp), parameter :: plot = 1 !Turn on/off saving a formatted file for plotting with e.g. python






end module parameters
