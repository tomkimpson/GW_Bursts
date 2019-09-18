module GravWaves

use parameters
use constants

implicit none


private FiniteDifference, Coeffs


public calc_GW 
       
      

contains



subroutine calc_GW(nrows,output)
!Arguments
integer(kind=dp) :: nrows
real(kind=dp), intent(in), dimension(nrows,entries+3) :: output

!Other
integer(kind=dp), parameter :: ncol = 7
real(kind=dp),dimension(nrows) :: m,r,theta,phi, vr,vtheta,vphi
real(kind=dp),dimension(nrows) :: t,x,y,z, vx,vy,vz
real(kind=dp), dimension(nrows,ncol) :: I, S1, S2, S3, M1, M2, M3
real(kind=dp), dimension(nrows,ncol-1) :: IDeriv, S1_Deriv,S2_Deriv, S3_Deriv, M1_Deriv, M2_Deriv,M3_Deriv, SDeriv, MDeriv, hbar
integer(kind=dp) :: order,j
real(kind=dp), parameter :: OBSTheta = 0.0_dp, OBSPhi = 0.0_dp
real(kind=dp),dimension(3) :: n
real(kind=dp), dimensiOn(nrows) :: hplus, hcross


!Load the data into a nice readable form
t = output(:,1)
r = output(:,2)
theta = output(:,3)
phi = output(:,4)
vr = output(:,13)
vtheta = output(:,14)
vphi = output(:,15)

!Covert to a cartesian form
m = sqrt(r**2 + a**2)
x = m*sin(theta)*cos(phi)
y = m*sin(theta)*sin(phi)
z = r*cos(theta)

vx = vr*(r*sin(theta)*cos(phi) / m) &
    +vtheta*(m*cos(theta)*cos(phi)) &
    - vphi*(sin(theta)*sin(phi))


vy = vr*(r*sin(theta)*sin(phi) / m) &
    +vtheta*(m*cos(theta)*sin(phi)) &
    + vphi*(sin(theta)*cos(phi))


vz = vr*(cos(theta)) &
    -vtheta*(sin(theta)) 


!Itensor
I(:,1) = t
I(:,2) = m0*x*x
I(:,3) = m0*y*y
I(:,4) = m0*z*z
I(:,5) = m0*x*y
I(:,6) = m0*x*z
I(:,7) = m0*y*z


do j = 1,nrows
!S tensor
S1(j,:) = vx(j) * I(j,:) 
S2(j,:) = vy(j) * I(j,:) 
S3(j,:) = vz(j) * I(j,:) 
!M tensor

M1(j,:) = x(j) * I(j,:) 
M2(j,:) = y(j) * I(j,:) 
M3(j,:) = z(j) * I(j,:) 
enddo

!Now calculate the derivatives via a finite difference scheme



!Second derivative of the inertia tensor
order = 2
call FiniteDifference(nrows,ncol,I,IDeriv,order)

!2nd derivatives of the S tensor
call FiniteDifference(nrows,ncol,S1,S1_Deriv,order)
call FiniteDifference(nrows,ncol,S2,S2_Deriv,order)
call FiniteDifference(nrows,ncol,S3,S3_Deriv,order)



!3rd derivative of M tensor
order = 3
call FiniteDifference(nrows,ncol,M1,M1_Deriv,order)
call FiniteDifference(nrows,ncol,M2,M2_Deriv,order)
call FiniteDifference(nrows,ncol,M3,M3_Deriv,order)


!Determine hplus and hcross

!Observer 'n' vector
n(1) = sin(OBSTheta)*cos(OBSPhi)
n(2) = sin(OBSTheta)*sin(OBSPhi)
n(3) = cos(OBSTheta)

SDeriv = n(1)*S1_Deriv + n(2)*S2_Deriv + n(3)*S3_Deriv 
MDeriv = n(1)*M1_Deriv + n(2)*M2_Deriv + n(3)*M3_Deriv 

hbar = 2.0_dp*(IDeriv - 2.0_dp*SDeriv + MDeriv)



  stop

end subroutine calc_GW





SUBROUTINE FiniteDifference(Nrows, Ncols,ArrayIn,ArrayOUT,m)

IMPLICIT NONE

integer(kind=dp) :: Nrows, Ncols !Length of array

real(kind=dp), dimension(Nrows,Ncols) :: ArrayIn !Define total input array
real(kind=dp), dimension(Nrows,Ncols-1) :: ArrayOut !Define total output array

integer(kind=dp),parameter :: nd = 5 !Number of grid points. Must be odd
integer(kind=dp),parameter :: nn = nd - 1 
integer(kind=dp) :: m

real(kind=dp), dimension(:), allocatable :: xdata
real(kind=dp), dimension(:,:), allocatable :: ydata

real(kind=dp) :: z
integer(kind=dp) :: i,j,k,l
real(kind=dp), dimension(Ncols-1) :: Iout


!Allocate temp arrays
ALLOCATE(xdata(0:nd-1))
ALLOCATE(ydata(0:nd-1,Ncols-1))


do j = 1,Nrows


!Get the target point
z = ArrayIn(j,1)

!Choose forward/backward/central modelling

if (j .LT. 1 + nn/2) then

!Forward modelling

do l=1,nd
xdata(l-1) = ArrayIn(j-1+l,1)
ydata(l-1,:) = ArrayIn(j-1+l,2:)
enddo



elseif (j + nn/2 .GT. Nrows) then
!Backward modelling


do l=1,nd
xdata(l-1) = ArrayIn(j-(nd-l),1)
ydata(l-1,:) = ArrayIn(j-(nd-l),2:)
enddo



else
!Central modelling

do l=1,nd
xdata(l-1) = ArrayIn(j-1-(nd/2)+l,1)
ydata(l-1,:) = ArrayIn(j-1-(nd/2)+l,2:)
enddo



endif

call Coeffs(nd,m,z,xdata,ydata,Iout)

ArrayOUT(j,:) = Iout


enddo
end subroutine FiniteDifference




subroutine Coeffs(nd,m,z,x,y,Iout)

implicit none
integer(kind=dp) :: nd, m,i,j,k
integer(kind=dp) :: nn,mn
real(kind=dp) :: z
real(kind=dp), dimension(0:nd-1) :: x
real(kind=dp), dimension(6) :: Iout
real(kind=dp), dimension(0:nd-1, 6) :: y
real(kind=dp), dimension(0:nd-1,0:m) :: c
real(kind=dp) :: c1,c2,c3,c4,c5



nn = nd-1
!Start
c1 = 1.0_dp
c4 = x(0) - z
c = 0.0_dp



c(0,0) = 1.0_dp

do i = 1,nn

mn = min(i,m)

c2 = 1.0_dp
c5 = c4
c4 = x(i) - z

do j = 0, i-1


    c3 = x(i) - x(j)
    c2 = c2*c3
    if ( j .EQ. i-1) then

        do k=mn,1,-1
            c(i,k) = c1*(k*c(i-1,k-1) - c5*c(i-1,k))/c2
        enddo
        c(i,0) = -c1*c5*c(i-1,0)/c2
    endif

    do k = mn,1,-1
        c(j,k) = (c4*c(j,k) - k*c(j,k-1))/c3
    enddo
    c(j,0) = c4*c(j,0)/c3
enddo

c1 = c2


enddo




do i = 1,6
Iout(i) =  SUM(c(:,m)*y(:,i))
enddo



end subroutine Coeffs
















end module GravWaves
