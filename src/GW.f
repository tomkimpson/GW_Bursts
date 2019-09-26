module GravWaves

use parameters
use constants

implicit none


private FiniteDifference, Coeffs


public calc_GW 
       
      

contains



subroutine calc_GW(nrows,output, hout)
!Arguments
integer(kind=dp) :: nrows
real(kind=dp), intent(in), dimension(nrows,entries+3) :: output
real(kind=dp),intent(out), dimension(nrows,2) :: hout

!Other
integer(kind=dp), parameter :: ncol = 7
real(kind=dp),dimension(nrows) :: m,r,theta,phi, vr,vtheta,vphi
real(kind=dp),dimension(nrows) :: t,x,y,z, vx,vy,vz
real(kind=dp), dimension(nrows,ncol) :: I, S1, S2, S3, M1, M2, M3
real(kind=dp), dimension(nrows,ncol-1) :: IDeriv, S1_Deriv,S2_Deriv, S3_Deriv, M1_Deriv, M2_Deriv,M3_Deriv, SDeriv, MDeriv,&
 hbar,h
real(kind=dp), dimension(nrows) :: hmag, hxx,hyy,hzz,hxy,hxz,hyz
integer(kind=dp) :: order,j,ll
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

vx = vr*r*sin(theta)*cos(phi) / m &
    +vtheta*m*cos(theta)*cos(phi) &
    - m*vphi*sin(theta)*sin(phi) 


vy = vr*r*sin(theta)*sin(phi) / m &
    +vtheta*m*cos(theta)*sin(phi) & 
    + m*vphi*sin(theta)*cos(phi) 


vz = vr*cos(theta) &
    -vtheta*r*sin(theta) 





!Itensor
I(:,1) = t
I(:,2) = m0*x*x
I(:,3) = m0*y*y
I(:,4) = m0*z*z
I(:,5) = m0*x*y
I(:,6) = m0*x*z
I(:,7) = m0*y*z











S1 = I
S2 = I
S3 = I
M1 = I
M2 = I
M3 = I

do j = 1,nrows
!S tensor
S1(j,2:ncol) = vx(j) * I(j,2:ncol) !not the 1st olumn which is time
S2(j,2:ncol) = vy(j) * I(j,2:ncol) 
S3(j,2:ncol) = vz(j) * I(j,2:ncol) 
!M tensor

M1(j,2:ncol) = x(j) * I(j,2:ncol) 
M2(j,2:ncol) = y(j) * I(j,2:ncol) 
M3(j,2:ncol) = z(j) * I(j,2:ncol) 
enddo

!Now calculate the derivatives via a finite difference scheme


!Second derivative of the inertia tensor
order = 2
call FiniteDifference(nrows,ncol,I,IDeriv,order)









!print *, S2(1,2),vy(1), I(1,2), vy(1)*I(1,2)






!print *, 'Checks', S2(1,:)
!stop




!2nd derivatives of the S tensor
call FiniteDifference(nrows,ncol,S1,S1_Deriv,order)
call FiniteDifference(nrows,ncol,S2,S2_Deriv,order)
call FiniteDifference(nrows,ncol,S3,S3_Deriv,order)








!3rd derivative of M tensor
order = 3
call FiniteDifference(nrows,ncol,M1,M1_Deriv,order)
call FiniteDifference(nrows,ncol,M2,M2_Deriv,order)
call FiniteDifference(nrows,ncol,M3,M3_Deriv,order)


!Observer 'n' vector
n(1) = sin(OBSTheta)*cos(OBSPhi)
n(2) = sin(OBSTheta)*sin(OBSPhi)
n(3) = cos(OBSTheta)

SDeriv = n(1)*S1_Deriv + n(2)*S2_Deriv + n(3)*S3_Deriv 
MDeriv = n(1)*M1_Deriv + n(2)*M2_Deriv + n(3)*M3_Deriv 

!Get hbar
hbar = 2.0_dp*(IDeriv - 2.0_dp*SDeriv + MDeriv)/OBSR

!print *, 'Ttot:', OBSR*hbar(1,:)/2.0_Dp
!stop






!do j = 1,5
!print*, IDeriv(j,1), SDeriv(j,1), MDeriv(j,1)
!enddo
!stop





!And convert to the non-trace reversed version
do j = 1,nrows
hmag(j) = hbar(j,1) + hbar(j,2) + hbar(j,3)

enddo


h = hbar

h(:,1) = h(:,1) - 0.50_dp*hmag(:)
h(:,2) = h(:,2) - 0.50_dp*hmag(:)
h(:,3) = h(:,3) - 0.50_dp*hmag(:)








!Andcalculate hplus andhcross, accounting for the half sign

hxx = h(:,1)
hyy = h(:,2)
hzz = h(:,3)
hxy = h(:,4)
hxz = h(:,5)
hyz = h(:,6)





hplus = 0.50_dp*(cos(OBSTheta)**2 * hxx - hyy + sin(OBSTheta)*hzz - sin(2.0_dp*OBSTheta)*hxz)
hcross = cos(OBSTheta)*hxy - sin(OBSTheta)*hyz




!And putput
hout(:,1) = hplus*OBSR/m0
hout(:,2) = hcross*OBSR/m0



!do j = 1,5
!print *, j, hplus(j), hplus(j)*OBSR/m0, OBSR/m0
!enddo






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
