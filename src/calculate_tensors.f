MODULE tensors

USE SHARED_CONSTANTS
IMPLICIT NONE

CONTAINS




subroutine ShrinkStepsize(errmax)
IMPLICIT NONE
real(kind=dp) :: htemp,errmax

htemp = S*h*errmax**Pshrink

h = sign(max(abs(htemp),0.10_dp*abs(h)),h)

end subroutine ShrinkStepsize


subroutine GrowStepsize(errmax)
IMPLICIT NONE
real(kind=dp) :: errmax, htemp


if (errmax .GT. errcon) then
htemp = S*h*errmax**Pgrow


else

htemp = 5.0_dp*h

endif


!if (htemp .LT. maxh) then
h = htemp
!endif





end subroutine GrowStepsize 






!SUBROUTINE FiniteDifference(NROWS,ArrayIN, ArrayOUT,m)


!Subroutine to calculate the m order derivatives of a time series
!Derivatives are made w.r.t the time series defined by ArrayIN

!implicit none

!integer(kind=dp) :: NROWS !Length of array
!real(kind=dp), dimension(NROWS) :: ArrayIN !Define IN array. A time 
!real(kind=dp), dimension(NROWS) :: ArrayOUT !Define OUT array


!integer(kind=dp) :: m ! This is the order of differentiation

!integer(kind=dp) :: j,l !Counting index
!real(kind=dp) :: z ! This is the point at which the derivative is to be evaluated
!integer(kind=dp),parameter :: nd = 5 !Number of grid points around the target point
!integer(kind=dp),parameter :: nn = nd - 1 


!real(kind=dp), dimension(0:nd-1) :: xdata !Temporary arrays just for storing grid points


!do j = 11,11

!Get the targe point
!z = ArrayIN(j) 

!print *, z

!Chose Forward/Backward/Central finite difference

!if (j .LT. 1+nn/2) then
!Forward

!do l=1,nd
!xdata(l-1) = ArrayIN(j-1+l)
!enddo



!elseif (j+nn/2 .GT. NROWS) then
!Backward

!do l = 1,nd
!xdata(l-1) = ArrayIN(j-(nd-l))
!enddo



!else
!central

!do l=1,nd
!xdata(l-1) = ArrayIN(j-1-(nd/2)+l)
!enddo








!print *, xdata




!endif





!print *, j


!enddo




!end subroutine FiniteDifference






!subroutine FiniteDifferenceCoeffs(nd,m,z)

!nd = number of grid points
!m = up to order. z = target point

!implicit none




!end subroutine FiniteDifferenceCoeffs










SUBROUTINE FiniteDifference(Nrows, Ncols,ArrayIn,ArrayOUT,m)

IMPLICIT NONE

integer(kind=dp) :: Nrows, Ncols !Length of array

real(kind=dp), dimension(Nrows,Ncols) :: ArrayIn !Define total input array
real(kind=dp), dimension(Nrows,Ncols-1) :: ArrayOut !Define total input array

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











subroutine GW2(N, Itensor,hOUT)

implicit none
integer(kind=dp) :: N !Length of array
real(kind=dp), dimension(N,6) :: ITensor
real(kind=dp), dimension(N,6) :: hii
real(kind=dp), dimension(N) :: hxx,hyy,hzz,hxy,hxz,hyz,hA, hB, hC
real(kind=dp), dimension(N) :: hxxB,hyyB,hzzB,hxyB,hxzB,hyzB,h

real(kind=dp), dimension(N,2) :: hOUT


hii = 2.0_dp*ITensor / OBSR


hxxB = hii(:,1)
hyyB = hii(:,2)
hzzB = hii(:,3)
hxyB = hii(:,4)
hxzB = hii(:,5)
hyzB = hii(:,6)


h = hxxB + hyyB + hzzB



hxx = hxxB - 0.50_dp*h
hyy = hyyB - 0.50_dp*h
hzz = hzzB - 0.50_dp*h
hxy = hxyB
hxz = hxzB
hyz = hyzB




print *, hxxB(1), hyyB(1), hzzB(1), hxxB(1) + hyyB(1) + hzzB(1)

!STOP







hA = Cos(OBSTheta)**2 * (hxx * Cos(OBSphi)**2 + hxy*Sin(2.0_dp*OBSphi) + hyy*Sin(OBSPhi)**2) + &
     hzz*Sin(OBSTheta)**2 - Sin(2.0_dp*OBSTheta)*(hxz *Cos(OBSPhi) + hyz*Sin(OBSPhi))


hB = Cos(OBSTheta)*(-0.50_dp*hxx*Sin(2.0_dp*OBSPhi)+hxy*Cos(2.0_dp*OBSPhi) + 0.50_dp*hyy*Sin(2.0_dp*OBSPhi)) + &
     Sin(OBSTheta)*(hxz*Sin(OBSPhi) - hyz*Cos(OBSPhi))



hC = hxx*Sin(OBSPhi)**2 - hxy*Sin(2.0_dp*OBSPhi) +hyy*Cos(OBSPhi)**2


!hOUT(:,1) = hA - hC
!hOUT(:,2) = 2.0_dp*hB



!print *, hA(1) - hC(1)


hOUT(:,1) = 0.50_dp*(cos(OBSTheta)**2 * hxx - hyy + sin(OBSTheta)*hzz - sin(2*OBSTheta)*hxz)
hOUT(:,2) = cos(OBSTheta)*hxy -sin(OBSTheta)*hyz










end subroutine GW2





SUBROUTINE GravWaves(r,theta,phi, &
                     rA,thetaA,phiA, &
                     rC,thetaC,phiC, &
                     hplus,hcross,&
                     out_array)

IMPLICIT NONE

real(kind=dp) :: r,theta,phi
real(kind=dp) :: rA,thetaA,phiA !First derivatives
real(kind=dp) :: rC,thetaC,phiC !Second derivatives
real(kind=dp) Ixx, Iyy,Izz, Ixy, Ixz, Iyz
real(kind=dp) :: hxx,hyy,hzz,hxy,hxz,hyz
real(kind=dp) :: hA, hB, hC
real(kind=dp) :: hplus, hcross
real(kind=dp), parameter :: m = m0
real(kind=dp), dimension(6) :: out_array
!These are the second derivatives of the mass quadrupole moment

Ixx =2.0_dp *m*(rA**2.0_dp *Cos(phi)**2.0_dp *Sin(theta)**2.0_dp  + r*Sin(theta)*(Cos(phi)**2.0_dp *(4.0_dp *rA*thetaA*Cos(theta) + rC*Sin(theta)) - 2.0_dp *phiA*rA*Sin(theta)*Sin(2.0_dp *phi)) + r**2.0_dp *(thetaA**2.0_dp *Cos(theta)**2.0_dp *Cos(phi)**2.0_dp  + thetaC*Cos(theta)*Cos(phi)**2.0_dp *Sin(theta) - (phiA**2.0_dp  + thetaA**2.0_dp )*Cos(phi)**2.0_dp *Sin(theta)**2.0_dp  - phiC*Cos(phi)*Sin(theta)**2.0_dp *Sin(phi) + phiA*(phiA*Sin(theta)**2.0_dp *Sin(phi)**2.0_dp  - thetaA*Sin(2.0_dp *theta)*Sin(2.0_dp *phi)))) 


Iyy =2.0_dp *m*(rA**2.0_dp *Sin(theta)**2.0_dp *Sin(phi)**2.0_dp  + r*Sin(phi)*(2.0_dp *rA*thetaA*Sin(2.0_dp *theta)*Sin(phi) + Sin(theta)**2.0_dp *(4.0_dp *phiA*rA*Cos(phi) + rC*Sin(phi))) + r**2.0_dp *(phiA**2.0_dp *Cos(phi)**2.0_dp *Sin(theta)**2.0_dp  + phiC*Cos(phi)*Sin(theta)**2.0_dp *Sin(phi) + thetaA**2.0_dp *Cos(theta)**2.0_dp *Sin(phi)**2.0_dp  + thetaC*Cos(theta)*Sin(theta)*Sin(phi)**2.0_dp  - phiA**2.0_dp *Sin(theta)**2.0_dp *Sin(phi)**2.0_dp  - thetaA**2.0_dp *Sin(theta)**2.0_dp *Sin(phi)**2.0_dp  + phiA*thetaA*Sin(2.0_dp *theta)*Sin(2.0_dp *phi)))

Izz =m*(2.0_dp *rA**2.0_dp *Cos(theta)**2.0_dp  + 2.0_dp *Cos(theta)*r*(rC*Cos(theta) - 4.0_dp *rA*thetaA*Sin(theta)) - r**2.0_dp *(2.0_dp *thetaA**2.0_dp *Cos(2.0_dp *theta) + thetaC*Sin(2.0_dp *theta)))

Ixy =m*(rA**2.0_dp *Sin(theta)**2.0_dp *Sin(2.0_dp *phi) + r*(4.0_dp *phiA*rA*Cos(phi)**2.0_dp *Sin(theta)**2.0_dp  + 2.0_dp *rA*thetaA*Sin(2.0_dp *theta)*Sin(2.0_dp *phi) + Sin(theta)**2.0_dp *(-4.0_dp *phiA*rA*Sin(phi)**2.0_dp  + rC*Sin(2.0_dp *phi))) + r**2.0_dp *(Cos(phi)**2.0_dp *Sin(theta)*(4.0_dp *phiA*thetaA*Cos(theta) + phiC*Sin(theta)) + thetaC*Cos(phi)*Sin(2.0_dp *theta)*Sin(phi) - Sin(theta)**2.0_dp *Sin(phi)*(2.0_dp *(2.0_dp *phiA**2.0_dp  + thetaA**2.0_dp )*Cos(phi) + phiC*Sin(phi)) + thetaA*(-2.0_dp *phiA*Sin(2.0_dp *theta)*Sin(phi)**2.0_dp  + thetaA*Cos(theta)**2.0_dp *Sin(2.0_dp *phi))))


Ixz =m*(rA**2.0_dp *Cos(phi)*Sin(2.0_dp *theta) + r*(4.0_dp *rA*thetaA*Cos(theta)**2.0_dp *Cos(phi) + Cos(phi)*(-4.0_dp *rA*thetaA*Sin(theta)**2.0_dp  + rC*Sin(2.0_dp *theta)) - 2.0_dp *phiA*rA*Sin(2.0_dp *theta)*Sin(phi)) + r**2.0_dp *(-(Cos(theta)*Sin(theta)*((phiA**2.0_dp  + 4.0_dp *thetaA**2.0_dp )*Cos(phi) + phiC*Sin(phi))) + Cos(theta)**2.0_dp *(thetaC*Cos(phi) - 2.0_dp *phiA*thetaA*Sin(phi)) + Sin(theta)**2.0_dp *(-(thetaC*Cos(phi)) + 2.0_dp *phiA*thetaA*Sin(phi))))

Iyz =m*(rA**2.0_dp *Sin(2.0_dp *theta)*Sin(phi) + r*(2.0_dp *phiA*rA*Cos(phi)*Sin(2.0_dp *theta) + (4.0_dp *rA*thetaA*Cos(2.0_dp *theta) + rC*Sin(2.0_dp *theta))*Sin(phi)) + r**2.0_dp *(Cos(theta)*Sin(theta)*(phiC*Cos(phi) - (phiA**2.0_dp  + 4.0_dp *thetaA**2.0_dp )*Sin(phi)) + Cos(theta)**2.0_dp *(2.0_dp *phiA*thetaA*Cos(phi) + thetaC*Sin(phi)) - Sin(theta)**2.0_dp *(2.0_dp *phiA*thetaA*Cos(phi) + thetaC*Sin(phi))))

 






!print *, 'Analytical Solution:', Ixx,Iyy,Izz,Ixy,Ixz,Iyz







!-----------------------------------------------

!Convert inertia tensor from natural units to SI
!Ixx = Msun*Ixx /(convert_m)**2
!Iyy = Msun*Iyy /(convert_m)**2
!Izz = Msun*Izz /(convert_m)**2
!Ixy = Msun*Ixy /(convert_m)**2
!Ixz = Msun*Ixz /(convert_m)**2
!Iyz = Msun*Iyz /(convert_m)**2





!print *, Ixx
!print *, OBSR
!print *, 2.0_dp*Newton_g/(OBSR*light_c**4) * Ixx
!print *, 'enddd here'
!STOP





!hxx = 2.0_dp*Newton_g*Ixx /(OBSR*light_c**4)
!hyy = 2.0_dp*Newton_g*Iyy /(OBSR*light_c**4)
!hzz = 2.0_dp*Newton_g*Izz /(OBSR*light_c**4)
!hxy = 2.0_dp*Newton_g*Ixy /(OBSR*light_c**4)
!hxz = 2.0_dp*Newton_g*Ixz /(OBSR*light_c**4)
!hyz = 2.0_dp*Newton_g*Iyz /(OBSR*light_c**4)


hxx = 2.0_dp/(OBSR) *Ixx
hyy = 2.0_dp/(OBSR) *Iyy
hzz = 2.0_dp/(OBSR) *Izz
hxy = 2.0_dp/(OBSR) *Ixy
hxz = 2.0_dp/(OBSR) *Ixz
hyz = 2.0_dp/(OBSR) *Iyz



!print *, Ixx, OBSR, hxx
!STOP


!
!-----------------------------------------------

!Now get the waveforms at different obs latitudes


!OBSTheta = 0.0_dp
call calc_h(hxx,hxy,hyy,hzz,hxz,hyz,hplus,hcross)
out_array(1) = hplus
out_array(2) = hcross
Printeger = 0
!OBSTheta = PI/2.0_dp

!OBSTheta = PI/4.0_dp
!call calc_h(hxx,hxy,hyy,hzz,hxz,hyz,hplus,hcross)
!out_array(3) = hplus
!out_array(4) = hcross

!OBSTheta = PI/2.0_dp
!call calc_h(hxx,hxy,hyy,hzz,hxz,hyz,hplus,hcross)
!out_array(5) = hplus
!out_array(6) = hcross



!----------------------------------


END SUBROUTINE GravWaves


subroutine calc_h(hxx,hxy,hyy,hzz,hxz,hyz,hplus,hcross)

IMPLICIT NONE

real(kind=dp) :: hxx,hxy,hyy,hzz,hxz,hyz
real(kind=dp) :: hA, hB, hC

real(kind=dp) :: hplus,hcross
real(kind=dp), parameter :: m = m0

hA = Cos(OBSTheta)**2 * (hxx * Cos(OBSphi)**2 + hxy*Sin(2.0_dp*OBSphi) + hyy*Sin(OBSPhi)**2) + &
     hzz*Sin(OBSTheta)**2 - Sin(2.0_dp*OBSTheta)*(hxz *Cos(OBSPhi) + hyz*Sin(OBSPhi))








hB = Cos(OBSTheta)*(-0.50_dp*hxx*Sin(2.0_dp*OBSPhi)+hxy*Cos(2.0_dp*OBSPhi) + 0.50_dp*hyy*Sin(2.0_dp*OBSPhi)) + &
     Sin(OBSTheta)*(hxz*Sin(OBSPhi) - hyz*Cos(OBSPhi))



hC = hxx*Sin(OBSPhi)**2 - hxy*Sin(2.0_dp*OBSPhi) +hyy*Cos(OBSPhi)**2



hplus = hA - hC
hcross = 2.0_dp*hB


!Normalise
hplus = hplus !* OBSR/m
hcross = hcross! * OBSR/m







end subroutine calc_h






SUBROUTINE calc_FourSpin(V0,V1,V2,V3,&
                         S0,S1,S2,S3,&
                         P0,P1,P2,P3,&
                         m0,&
                         Sdot_0, Sdot_1, Sdot_2, Sdot_3)

IMPLICIT NONE
real (kind=dp) :: V0,V1,V2,V3
real(kind=dp) :: S0,S1,S2,S3 
REAL(kind=dp) :: P0,P1,P2,P3
real(kind=dp) :: Scoeff, m0
real(kind=dp) :: Sdot2_0, Sdot2_1, Sdot2_2, Sdot2_3
real(kind=dp) :: Sdot1_0, Sdot1_1, Sdot1_2, Sdot1_3                                         
real(kind=dp) :: Sdot_0, Sdot_1, Sdot_2, Sdot_3

Scoeff = -(R_0212*S0*V2*S_12-R_1213*S2*V1*S_13-R_1212*S2*V1*S_12-R_1201*S2*V1*S_01-R_2302   &
 *S3*V2*S_02+R_0312*S0*V3*S_12+R_0313*S0*V3*S_13+R_0302*S0*V3*S_02-R_0323*S3*V0*S_23+R_2301  &
 *S2*V3*S_01-R_1303*S3*V1*S_03-R_1202*S2*V1*S_02-R_1203*S2*V1*S_03-R_1223*S2*V1*S_23-R_0301  &
 *S3*V0*S_01-R_0312*S3*V0*S_12-R_1302*S3*V1*S_02-R_1301*S3*V1*S_01+R_2323*S2*V3*S_23-R_0303  &
 *S3*V0*S_03-R_2323*S3*V2*S_23-R_0313*S3*V0*S_13-R_0302*S3*V0*S_02-R_2312*S3*V2*S_12-R_2301  &
 *S3*V2*S_01-R_2313*S3*V2*S_13+R_0113*S0*V1*S_13+R_0213*S0*V2*S_13+R_0223*S0*V2*S_23+R_0112  &
 *S0*V1*S_12-R_0102*S1*V0*S_02+R_0203*S0*V2*S_03-R_0112*S1*V0*S_12-R_0123*S1*V0*S_23+R_0103  &
 *S0*V1*S_03+R_0102*S0*V1*S_02-R_1313*S3*V1*S_13-R_1312*S3*V1*S_12-R_1323*S3*V1*S_23-R_2303  &
 *S3*V2*S_03-R_0113*S1*V0*S_13-R_0101*S1*V0*S_01-R_0103*S1*V0*S_03+R_1223*S1*V2*S_23+R_1213  &
 *S1*V2*S_13+R_1302*S1*V3*S_02+R_0303*S0*V3*S_03+R_0301*S0*V3*S_01+R_0323*S0*V3*S_23-R_0223  &
 *S2*V0*S_23+R_1203*S1*V2*S_03+R_1202*S1*V2*S_02+R_1201*S1*V2*S_01+R_1323*S1*V3*S_23+R_1303  &
 *S1*V3*S_03+R_1312*S1*V3*S_12+R_1212*S1*V2*S_12-R_0213*S2*V0*S_13-R_0201*S2*V0*S_01-R_0202  &
 *S2*V0*S_02-R_0212*S2*V0*S_12+R_1313*S1*V3*S_13+R_1301*S1*V3*S_01-R_0203*S2*V0*S_03+R_2303  &
 *S2*V3*S_03+R_0101*S0*V1*S_01+R_2302*S2*V3*S_02+R_2313*S2*V3*S_13+R_0123*S0*V1*S_23+R_2312  &
 *S2*V3*S_12+R_0202*S0*V2*S_02+R_0201*S0*V2*S_01)/m0**2




Sdot2_0 = Scoeff*P0
Sdot2_1 = Scoeff*P1
Sdot2_2 = Scoeff*P2
Sdot2_3 = Scoeff*P3


Sdot1_0 = -G0_00*V0*S0-G0_01*V0*S1-G0_02*V0*S2-G0_03*V0*S3-G0_01*V1*S0-G0_11*V1*S1-G0_12*V1 &
 *S2-G0_13*V1*S3-G0_02*V2*S0-G0_12*V2*S1-G0_22*V2*S2-G0_23*V2*S3-G0_03*V3*S0-G0_13*V3*S1     &
-G0_23*V3*S2-G0_33*V3*S3

Sdot1_1 = -G1_00*V0*S0-G1_01*V0*S1-G1_02*V0*S2-G1_03*V0*S3-G1_01*V1*S0-G1_11*V1*S1-G1_12*V1 &
 *S2-G1_13*V1*S3-G1_02*V2*S0-G1_12*V2*S1-G1_22*V2*S2-G1_23*V2*S3-G1_03*V3*S0-G1_13*V3*S1     &
-G1_23*V3*S2-G1_33*V3*S3

Sdot1_2 = -G2_00*V0*S0-G2_01*V0*S1-G2_02*V0*S2-G2_03*V0*S3-G2_01*V1*S0-G2_11*V1*S1-G2_12*V1 &
 *S2-G2_13*V1*S3-G2_02*V2*S0-G2_12*V2*S1-G2_22*V2*S2-G2_23*V2*S3-G2_03*V3*S0-G2_13*V3*S1     &
-G2_23*V3*S2-G2_33*V3*S3

Sdot1_3 = -G3_00*V0*S0-G3_01*V0*S1-G3_02*V0*S2-G3_03*V0*S3-G3_01*V1*S0-G3_11*V1*S1-G3_12*V1 &
 *S2-G3_13*V1*S3-G3_02*V2*S0-G3_12*V2*S1-G3_22*V2*S2-G3_23*V2*S3-G3_03*V3*S0-G3_13*V3*S1     &
-G3_23*V3*S2-G3_33*V3*S3



Sdot_0 = Sdot1_0 + lambda*Sdot2_0
Sdot_1 = Sdot1_1 + lambda*Sdot2_1
Sdot_2 = Sdot1_2 + lambda*Sdot2_2
Sdot_3 = Sdot1_3 + lambda*Sdot2_3





END SUBROUTINE calc_FourSpin




SUBROUTINE calc_FourMom(V0,V1,V2,V3, &
                        P0,P1,P2,P3, &
                        Pdot_0, Pdot_1, Pdot_2, Pdot_3)
implicit none
real (kind=dp) :: V0,V1,V2,V3
real(kind=dp) :: P0,P1,P2,P3 
real(kind=dp) :: Pdot1_0, Pdot1_1, Pdot1_2, Pdot1_3
real(kind=dp) :: RS0_0, RS0_1, RS0_2, RS0_3, RS1_0, RS1_1, RS1_2, RS1_3,                     &
                 RS2_0, RS2_1, RS2_2, RS2_3, RS3_0, RS3_1, RS3_2, RS3_3
real(kind=dp) :: Pdot2_0, Pdot2_1, Pdot2_2, Pdot2_3                                         
real(kind=dp) :: Pdot_0, Pdot_1, Pdot_2, Pdot_3



Pdot1_0 = -G0_00*V0*P0-G0_01*V0*P1-G0_02*V0*P2-G0_03*V0*P3-G0_01*V1*P0-G0_11*V1*P1-G0_12*V1 &
 *P2-G0_13*V1*P3-G0_02*V2*P0-G0_12*V2*P1-G0_22*V2*P2-G0_23*V2*P3-G0_03*V3*P0-G0_13*V3*P1     &
-G0_23*V3*P2-G0_33*V3*P3

Pdot1_1 = -G1_00*V0*P0-G1_01*V0*P1-G1_02*V0*P2-G1_03*V0*P3-G1_01*V1*P0-G1_11*V1*P1-G1_12*V1 &
 *P2-G1_13*V1*P3-G1_02*V2*P0-G1_12*V2*P1-G1_22*V2*P2-G1_23*V2*P3-G1_03*V3*P0-G1_13*V3*P1     &
-G1_23*V3*P2-G1_33*V3*P3

Pdot1_2 = -G2_00*V0*P0-G2_01*V0*P1-G2_02*V0*P2-G2_03*V0*P3-G2_01*V1*P0-G2_11*V1*P1-G2_12*V1 &
 *P2-G2_13*V1*P3-G2_02*V2*P0-G2_12*V2*P1-G2_22*V2*P2-G2_23*V2*P3-G2_03*V3*P0-G2_13*V3*P1     &
-G2_23*V3*P2-G2_33*V3*P3

Pdot1_3 = -G3_00*V0*P0-G3_01*V0*P1-G3_02*V0*P2-G3_03*V0*P3-G3_01*V1*P0-G3_11*V1*P1-G3_12*V1 &
 *P2-G3_13*V1*P3-G3_02*V2*P0-G3_12*V2*P1-G3_22*V2*P2-G3_23*V2*P3-G3_03*V3*P0-G3_13*V3*P1     &
-G3_23*V3*P2-G3_33*V3*P3





RS0_0 = R0_001*S_01 + R0_002*S_02 + R0_003*S_03 + R0_012*S_12 + R0_013*S_13 + R0_023*S_23
RS0_1 = R0_101*S_01 + R0_102*S_02 + R0_103*S_03 + R0_112*S_12 + R0_113*S_13 + R0_123*S_23
RS0_2 = R0_201*S_01 + R0_202*S_02 + R0_203*S_03 + R0_212*S_12 + R0_213*S_13 + R0_223*S_23
RS0_3 = R0_301*S_01 + R0_302*S_02 + R0_303*S_03 + R0_312*S_12 + R0_313*S_13 + R0_323*S_23

RS1_0 = R1_001*S_01 + R1_002*S_02 + R1_003*S_03 + R1_012*S_12 + R1_013*S_13 + R1_023*S_23
RS1_1 = R1_101*S_01 + R1_102*S_02 + R1_103*S_03 + R1_112*S_12 + R1_113*S_13 + R1_123*S_23
RS1_2 = R1_201*S_01 + R1_202*S_02 + R1_203*S_03 + R1_212*S_12 + R1_213*S_13 + R1_223*S_23
RS1_3 = R1_301*S_01 + R1_302*S_02 + R1_303*S_03 + R1_312*S_12 + R1_313*S_13 + R1_323*S_23

RS2_0 = R2_001*S_01 + R2_002*S_02 + R2_003*S_03 + R2_012*S_12 + R2_013*S_13 + R2_023*S_23
RS2_1 = R2_101*S_01 + R2_102*S_02 + R2_103*S_03 + R2_112*S_12 + R2_113*S_13 + R2_123*S_23
RS2_2 = R2_201*S_01 + R2_202*S_02 + R2_203*S_03 + R2_212*S_12 + R2_213*S_13 + R2_223*S_23
RS2_3 = R2_301*S_01 + R2_302*S_02 + R2_303*S_03 + R2_312*S_12 + R2_313*S_13 + R2_323*S_23

RS3_0 = R3_001*S_01 + R3_002*S_02 + R3_003*S_03 + R3_012*S_12 + R3_013*S_13 + R3_023*S_23
RS3_1 = R3_101*S_01 + R3_102*S_02 + R3_103*S_03 + R3_112*S_12 + R3_113*S_13 + R3_123*S_23
RS3_2 = R3_201*S_01 + R3_202*S_02 + R3_203*S_03 + R3_212*S_12 + R3_213*S_13 + R3_223*S_23
RS3_3 = R3_301*S_01 + R3_302*S_02 + R3_303*S_03 + R3_312*S_12 + R3_313*S_13 + R3_323*S_23



Pdot2_0 = -(RS0_0*V0 + RS0_1*V1 + RS0_2*V2 + RS0_3*V3)

Pdot2_1 = -(RS1_0*V0 + RS1_1*V1 + RS1_2*V2 + RS1_3*V3)

Pdot2_2 = -(RS2_0*V0 + RS2_1*V1 + RS2_2*V2 + RS2_3*V3)

Pdot2_3 = -(RS3_0*V0 + RS3_1*V1 + RS3_2*V2 + RS3_3*V3)


Pdot_0 = Pdot1_0 + lambda*Pdot2_0
Pdot_1 = Pdot1_1 + lambda*Pdot2_1
Pdot_2 = Pdot1_2 + lambda*Pdot2_2
Pdot_3 = Pdot1_3 + lambda*Pdot2_3

END SUBROUTINE calc_FourMom









SUBROUTINE calc_FourVel(m0,P0,P1,P2,P3,&
                           V0,V1,V2,V3,&
                           g00,g11,g22,g33,g30)
IMPLICIT NONE
real(kind = dp) m0, P0,P1,P2,P3, M
real(kind=dp) V0,V1,V2,V3, PV
real(kind=dp) Vsq,g00,g11,g22,g33,g30
M = 1.00_dp

Delta = 2.0_dp*(m0**2 + lambda*                                                             &
(R_0101*S_01**2+R_0202*S_02**2+R_2323*S_23**2+R_0203*S_02*S_03+R_1301*S_13*S_01             & 
+R_1201*S_12*S_01+R_0313*S_03*S_13+R_0312*S_03*S_12+R_0113*S_01*S_13+R_0112*S_01*S_12       &
+R_0301*S_03*S_01+R_1212*S_12**2+R_0201*S_02*S_01+R_0303*S_03**2+R_0302*S_03*S_02           &
+R_1312*S_13*S_12+R_0103*S_01*S_03+R_0102*S_01*S_02+R_1313*S_13**2+R_0123*S_01*S_23         &
+R_1202*S_12*S_02+R_1203*S_12*S_03+R_2313*S_23*S_13+R_2312*S_23*S_12+R_1302*S_13*S_02       &
+R_2303*S_23*S_03+R_2302*S_23*S_02+R_0213*S_02*S_13+R_0212*S_02*S_12+R_2301*S_23*S_01       &
+R_1303*S_13*S_03+R_1213*S_12*S_13+R_1223*S_12*S_23+R_0223*S_02*S_23+R_1323*S_13*S_23       &
+R_0323*S_03*S_23)/4.0_dp)


RPS_0 = S_01**2*R_1301*P3+S_01**2*R_1201*P2+S_02*R_2312*P3*S_12+S_01*R_1203*P2*S_03         &
+S_02*R_2313*P3*S_13-S_01**2*R_0101*P0+S_01*R_1213*P2*S_13-S_02*R_0203*P0*S_03-S_02*R_1201  &
 *P1*S_01-S_03*R_2301*P2*S_01-S_03*R_1301*P1*S_01-S_03*R_0301*P0*S_01-S_02*R_0201*P0*S_01    &
-S_02*R_1203*P1*S_03-S_01*R_0103*P0*S_03-S_03*R_2302*P2*S_02-S_03*R_1302*P1*S_02-S_01       &
 *R_0112*P0*S_12-S_03*R_2312*P2*S_12-S_03*R_1312*P1*S_12-S_03*R_0312*P0*S_12-S_02*R_1212*P1  &
 *S_12-S_02*R_0212*P0*S_12-S_01*R_0102*P0*S_02-S_01*R_0113*P0*S_13-S_03*R_2313*P2*S_13       &
-S_03*R_1313*P1*S_13-S_03*R_0313*P0*S_13-S_02*R_1213*P1*S_13-S_02*R_0213*P0*S_13-S_02       &
 *R_1223*P1*S_23-S_01*R_0123*P0*S_23-S_03*R_2323*P2*S_23-S_03*R_1323*P1*S_23+S_01*R_1302*P3  &
 *S_02+S_01*R_1202*P2*S_02+S_01*R_1313*P3*S_13+S_01*R_1312*P3*S_12+S_01*R_1212*P2*S_12       &
+S_02**2*R_2302*P3-S_03*R_0302*P0*S_02-S_02**2*R_0202*P0-S_02**2*R_1202*P1-S_03**2*R_0303   &
 *P0-S_03**2*R_1303*P1-S_03**2*R_2303*P2+S_02*R_2323*P3*S_23-S_03*R_0323*P0*S_23+S_01*R_1323 &
 *P3*S_23+S_01*R_1223*P2*S_23-S_02*R_0223*P0*S_23+S_02*R_2301*P3*S_01+S_01*R_1303*P3*S_03    &
+S_02*R_2303*P3*S_03

RPS_1 = S_12**2*R_2312*P3+S_12*R_2303*P3*S_03+S_12*R_2313*P3*S_13+S_12*R_2323*P3*S_23       &
-S_01**2*R_0101*P1-S_01**2*R_0201*P2-S_01**2*R_0301*P3-S_12**2*R_0212*P0-S_13**2*R_1313*P1  &
-S_13**2*R_0313*P0-S_01*R_0103*P1*S_03-S_13*R_2301*P2*S_01-S_13*R_0301*P0*S_01-S_13*R_1312  &
 *P1*S_12-S_13*R_0312*P0*S_12-S_13*R_2312*P2*S_12-S_12*R_0201*P0*S_01-S_13*R_2302*P2*S_02    &
-S_13*R_1302*P1*S_02-S_12*R_1203*P1*S_03-S_12*R_0203*P0*S_03-S_12*R_1201*P1*S_01-S_01       &
 *R_0303*P3*S_03-S_01*R_0203*P2*S_03-S_12**2*R_1212*P1-S_01*R_0323*P3*S_23-S_01*R_0223*P2    &
 *S_23-S_01*R_0123*P1*S_23-S_12*R_1223*P1*S_23-S_12*R_0223*P0*S_23-S_13*R_0323*P0*S_23       &
-S_13*R_2323*P2*S_23-S_13*R_1323*P1*S_23+S_12*R_2302*P3*S_02-S_12*R_1202*P1*S_02-S_01       &
 *R_0202*P2*S_02-S_13*R_0302*P0*S_02-S_01*R_0113*P1*S_13-S_01*R_0302*P3*S_02-S_12*R_0202     &
 *P0*S_02-S_01*R_0102*P1*S_02-S_12*R_0213*P0*S_13-S_01*R_0313*P3*S_13-S_01*R_0213*P2*S_13    &
-S_13**2*R_2313*P2-S_01*R_0212*P2*S_12-S_01*R_0112*P1*S_12-S_12*R_1213*P1*S_13-S_01*R_0312  &
 *P3*S_12-S_13*R_1303*P1*S_03-S_13*R_0303*P0*S_03-S_13*R_1301*P1*S_01-S_13*R_2303*P2*S_03    &
+S_12*R_2301*P3*S_01

RPS_2 = -S_12*R_1201*P2*S_01-S_02*R_0301*P3*S_01-S_02*R_0201*P2*S_01-S_23*R_1312*P1*S_12    &
-S_02*R_0101*P1*S_01-S_02*R_0112*P1*S_12-S_02*R_0312*P3*S_12+S_12**2*R_0112*P0-S_12*R_1203  &
 *P2*S_03+S_12*R_0101*P0*S_01-S_23*R_0312*P0*S_12+S_12*R_0103*P0*S_03-S_23*R_2312*P2*S_12    &
+S_12*R_0113*P0*S_13-S_02**2*R_0102*P1-S_12*R_1302*P3*S_02-S_23*R_1303*P1*S_03+S_12*R_0102  &
 *P0*S_02-S_02**2*R_0302*P3-S_23*R_0303*P0*S_03-S_12*R_1313*P3*S_13-S_12**2*R_1312*P3        &
-S_12**2*R_1212*P2-S_23**2*R_0323*P0-S_23*R_2303*P2*S_03-S_23**2*R_1323*P1-S_02**2*R_0202   &
 *P2-S_23*R_0301*P0*S_01-S_12*R_1301*P3*S_01-S_23*R_2301*P2*S_01-S_23*R_1301*P1*S_01-S_02    &
 *R_0212*P2*S_12-S_12*R_1303*P3*S_03-S_02*R_0323*P3*S_23-S_23**2*R_2323*P2-S_02*R_0303*P3    &
 *S_03+S_12*R_0123*P0*S_23-S_12*R_1223*P2*S_23-S_12*R_1213*P2*S_13-S_02*R_0203*P2*S_03-S_23  &
 *R_2302*P2*S_02-S_02*R_0103*P1*S_03-S_23*R_1302*P1*S_02-S_23*R_0302*P0*S_02-S_23*R_2313*P2  &
 *S_13-S_23*R_1313*P1*S_13-S_02*R_0113*P1*S_13-S_02*R_0313*P3*S_13-S_02*R_0213*P2*S_13       &
-S_12*R_1202*P2*S_02-S_12*R_1323*P3*S_23-S_23*R_0313*P0*S_13-S_02*R_0223*P2*S_23-S_02       &
 *R_0123*P1*S_23

RPS_3 = S_23*R_0201*P0*S_01-S_13*R_1301*P3*S_01-S_03*R_0201*P2*S_01-S_03*R_0101*P1*S_01     &
-S_13*R_1201*P2*S_01-S_03*R_0301*P3*S_01-S_23*R_2303*P3*S_03-S_13*R_1303*P3*S_03+S_23       &
 *R_1212*P1*S_12+S_13*R_0103*P0*S_03+S_23*R_0212*P0*S_12-S_23*R_2312*P3*S_12-S_13*R_1203*P2  &
 *S_03+S_13*R_0123*P0*S_23-S_13*R_1323*P3*S_23-S_13*R_1223*P2*S_23-S_03*R_0323*P3*S_23       &
-S_03*R_0123*P1*S_23-S_03*R_0223*P2*S_23+S_23**2*R_0223*P0+S_23**2*R_1223*P1-S_23*R_2302*P3 &
 *S_02+S_13*R_0102*P0*S_02+S_23*R_1202*P1*S_02-S_23**2*R_2323*P3+S_23*R_0202*P0*S_02-S_13    &
 *R_1202*P2*S_02-S_13*R_1302*P3*S_02-S_03*R_0313*P3*S_13+S_13**2*R_0113*P0-S_03*R_0213*P2    &
 *S_13-S_03*R_0102*P1*S_02+S_23*R_0213*P0*S_13-S_13*R_1312*P3*S_12-S_13*R_1212*P2*S_12-S_03  &
 *R_0302*P3*S_02-S_03*R_0202*P2*S_02-S_03*R_0113*P1*S_13+S_23*R_1213*P1*S_13-S_13**2*R_1213  &
 *P2-S_03*R_0212*P2*S_12-S_03*R_0112*P1*S_12-S_23*R_2313*P3*S_13+S_13*R_0112*P0*S_12-S_13**2 &
 *R_1313*P3-S_03*R_0312*P3*S_12+S_23*R_1203*P1*S_03+S_23*R_0203*P0*S_03+S_13*R_0101*P0*S_01  &
+S_23*R_1201*P1*S_01-S_03**2*R_0303*P3-S_23*R_2301*P3*S_01-S_03**2*R_0203*P2-S_03**2*R_0103 &
 *P1




!CHOOSE PARAMETERIZATION. Use this if you want the parameterizationisproper time

V0 = - (P0 + lambda*RPS_0/Delta)/m0**2
V1 = - (P1 + lambda*RPS_1/Delta)/m0**2
V2 = - (P2 + lambda*RPS_2/Delta)/m0**2
V3 = - (P3 + lambda*RPS_3/Delta)/m0**2


Vsq = g00*V0**2 + g11*V1**2 + g22*V2**2 + g33*V3**2 + 2.0_dp*g30*V0*V3


PV = -sqrt(-1.0_dp/Vsq)

V0 = PV*V0
V1 = PV*V1
V2 = PV*V2
V3 = PV*V3



!This is a check that = -1 if the parameterization is the proper time
!print *,'CHECK:', g00*V0**2 + g11*V1**2 + g22*V2**2 + g33*V3**2 + 2.0_dp*g30*V0*V3
!
!STOP

!-----NEW------!



!--------------ORIG----------!

!V0 = 1.0_dp

!PV = -M**2*V0/(P0 + lambda*RPS_0/Delta)

!V1 = -PV/M**2*(P1 + lambda*RPS_1/Delta)

!V2 = -PV/M**2*(P2 + lambda*RPS_2/Delta)


!print *,'Orig V2',V2, lambda, RPS_2,Delta

!V3 = -PV/M**2*(P3 + lambda*RPS_3/Delta)


!--------------ORIG----------!


END SUBROUTINE calc_FourVel




!SUBROUTINE calc_FourAccel(m0,V0,V1,V2,V3)
!IMPLICIT NONE
!real(kind = dp) m0, M
!real(kind=dp) a0,a1,a2,a3
!real(kind=dp) V0,V1,V2,V3

!real(kind=dp) :: a00,a01,a02,a03



!M = 1.00_dp

!a00 = (&
 !1     R_0000*V0*S_00 + &
 !     R_0001*V0*S_01 + &
 !     R_0002*V0*S_02 + &
 !     R_0003*V0*S_03 + &
 !     R_0010*V0*S_10 + &
 !     R_0011*V0*S_11 + &
!      R_0012*V0*S_12 + &
!      R_0013*V0*S_13 + &
!      R_0020*V0*S_20 + &
!      R_0021*V0*S_21 + &
!      R_0022*V0*S_22 + &
!      R_0023*V0*S_23 + &
!      R_0030*V0*S_30 + &
   !   R_0031*V0*S_31 + &
  !    R_0032*V0*S_32 + &
 !     R_0033*V0*S_33 )


!a01 = (&
     ! R_0100*V1*S_00 + &
     ! R_0101*V1*S_01 + &
     ! R_0102*V1*S_02 + &
     ! R_0103*V1*S_03 + &
     ! R_0110*V1*S_10 + &
     ! R_0111*V1*S_11 + &
     ! R_0112*V1*S_12 + &
     ! R_0113*V1*S_13 + &
     ! R_0120*V1*S_20 + &
     ! R_0121*V1*S_21 + &
     ! R_0122*V1*S_22 + &
     ! R_0123*V1*S_23 + &
     ! R_0130*V1*S_30 + &
 !     R_0131*V1*S_31 + &
 !     R_0132*V1*S_32 + &
 !     R_0133*V1*S_33 )


!a02 = (&
!      R_0200*V2*S_00 + &
!      R_0201*V2*S_01 + &
 !     R_0202*V2*S_02 + &
!      R_0203*V2*S_03 + &
!      R_0210*V2*S_10 + &
!      R_0211*V2*S_11 + &
!      R_0212*V2*S_12 + &
 !     R_0213*V2*S_13 + &
!      R_0220*V2*S_20 + &
!      R_0221*V2*S_21 + &
!      R_0222*V2*S_22 + &
!      R_0223*V2*S_23 + &
!      R_0230*V2*S_30 + &
   !   R_0231*V2*S_31 + &
  !    R_0232*V2*S_32 + &
 !     R_0233*V2*S_33 )

!  end SUBROUTINE calc_FourAccel



SUBROUTINE calc_christoffel(r,theta)
real(kind=dp) :: r,theta, Sg, Dl, M


M=1.00_dp
Sg = r**2+a**2*cos(theta)**2  
Dl = r**2 + a**2 - 2.0_dp*M*r

G0_00 = 0.0_dp
G0_01 = M*(r**2+a**2)*(r**2-a**2*cos(theta)**2)/(Sg**2*Dl)
G0_02 = -2.0_dp*M*r*a**2*sin(theta)*cos(theta)/Sg**2
G0_03 = 0.0_dp
G0_11 = 0.0_dp
G0_12 = 0.0_dp
G0_13 = -M*a*((3.0_dp*r**2-a**2)*(r**2+a**2)-a**2*(r**2-a**2)*sin(theta)**2)        &
*sin(theta)**2/(Sg**2*Dl)
G0_22 = 0.0_dp
G0_23 = 2.0_dp*M*r*a**3*sin(theta)**3*cos(theta)/Sg**2
G0_33 = 0.0_dp

G1_00 = M*Dl*(r**2-a**2*cos(theta)**2)/Sg**3
G1_01 = 0.0_dp
G1_02 = 0.0_dp
G1_03 = -M*a*Dl*(r**2-a**2*cos(theta)**2)*sin(theta)**2/Sg**3
G1_11 = (-M*(r**2-a**2)+a**2*sin(theta)**2*(r-M))/(Sg*Dl)
G1_12 = -a**2*sin(theta)*cos(theta)/Sg
G1_13 = 0.0_dp
G1_22 = -r*Dl/Sg
G1_23 = 0.0_dp
G1_33 = -(r*(a**2+r**2)**2-a**2*sin(theta)**2*(-a**2*sin(theta)**2*(r-M)            &
-M*a**2+r*(2.0_dp*(r**2+a**2)+M*r)))*Dl*sin(theta)**2/Sg**3

G2_00 = -2.0_dp*M*r*a**2*sin(theta)*cos(theta)/Sg**3
G2_01 = 0.0_dp
G2_02 = 0.0_dp
G2_03 = 2.0_dp*M*r*a*(r**2+a**2)*sin(theta)*cos(theta)/Sg**3
G2_11 = a**2*sin(theta)*cos(theta)/(Sg*Dl)
G2_12 = r/Sg
G2_13 = 0.0_dp
G2_22 = -a**2*sin(theta)*cos(theta)/Sg
G2_23 = 0.0_dp
G2_33 = -((r**2+a**2)**3-a**2*Dl*(2.0_dp*(r**2+a**2)-a**2*sin(theta)**2)*sin(theta)**2) &
*sin(theta)*cos(theta)/Sg**3

G3_00 = 0.0_dp
G3_01 = M*a*(r**2-a**2*cos(theta)**2)/(Sg**2*Dl)
G3_02 = -2.0_dp*M*r*a*(cos(theta)/sin(theta))/Sg**2
G3_03 = 0.0_dp
G3_11 = 0.0_dp
G3_12 = 0.0_dp
G3_13 = (r*Dl*(r**2+a**2)+a**4*(r-M)*sin(theta)**4-a**2*(a**2+r**2)*(2.0_dp*r-M)*sin(theta)**2) &
/(Sg**2*Dl)
G3_22 = 0.0_dp
G3_23 = ((r**2+a**2)**2+sin(theta)**4*a**4-2.0_dp*a**2*(r**2-M*r+a**2)*sin(theta)**2)*(cos(theta)/sin(theta))/Sg**2 
G3_33 = 0.0_dp





END SUBROUTINE calc_christoffel





SUBROUTINE calc_spin_antisym(r,theta, &
                             S0,S1,S2,S3,&
                             P0,P1,P2,P3,&
                             H00,H01,H02,H03,&
                             H10,H11,H12,H13,&
                             H20,H21,H22,H23,&
                             H30,H31, H32,H33)

IMPLICIT NONE
real(kind=dp) :: r,theta
real(kind=dp) ::  S0,S1,S2,S3,&
                  P0,P1,P2,P3,&
                  H00,H01,H02,H03,&
                  H10,H11,H12,H13,&
                  H20,H21,H22,H23,&
                  H30,H31, H32,H33

real(kind=dp) :: Sg, eta0123


Sg = r**2+a**2*cos(theta)**2  
eta0123 = Sg*sin(theta)


S_01 = eta0123/m0*((H00*H11-H01*H01)*(P2*S3-P3*S2)                    &
+ (H00*H12-H01*H02)*(P3*S1-P1*S3) + (H00*H13-H01*H03)*(P1*S2-P2*S1)   & 
+ (H10*H12-H11*H02)*(P0*S3-P3*S0) + (H10*H13-H11*H03)*(P2*S0-P0*S2)   &
+ (H20*H13-H21*H03)*(P0*S1-P1*S0))

S_02 = eta0123/m0*((H00*H21-H02*H01)*(P2*S3-P3*S2)                    &
+ (H00*H22-H02*H02)*(P3*S1-P1*S3) + (H00*H23-H02*H03)*(P1*S2-P2*S1)   & 
+ (H10*H22-H12*H02)*(P0*S3-P3*S0) + (H10*H23-H12*H03)*(P2*S0-P0*S2)   &
+ (H20*H23-H22*H03)*(P0*S1-P1*S0))

S_03 = eta0123/m0*((H00*H31-H03*H01)*(P2*S3-P3*S2)                    &
+ (H00*H32-H03*H02)*(P3*S1-P1*S3) + (H00*H33-H03*H03)*(P1*S2-P2*S1)   & 
+ (H10*H32-H13*H02)*(P0*S3-P3*S0) + (H10*H33-H13*H03)*(P2*S0-P0*S2)   &
+ (H20*H33-H23*H03)*(P0*S1-P1*S0))

S_12 = eta0123/m0*((H01*H21-H02*H11)*(P2*S3-P3*S2)                    &
+ (H01*H22-H02*H12)*(P3*S1-P1*S3) + (H01*H23-H02*H13)*(P1*S2-P2*S1)   & 
+ (H11*H22-H12*H12)*(P0*S3-P3*S0) + (H11*H23-H12*H13)*(P2*S0-P0*S2)   & 
+ (H21*H23-H22*H13)*(P0*S1-P1*S0))

S_13 = eta0123/m0*((H01*H31-H03*H11)*(P2*S3-P3*S2)                    &
+ (H01*H32-H03*H12)*(P3*S1-P1*S3) + (H01*H33-H03*H13)*(P1*S2-P2*S1)   & 
+ (H11*H32-H13*H12)*(P0*S3-P3*S0) + (H11*H33-H13*H13)*(P2*S0-P0*S2)   &
+ (H21*H33-H23*H13)*(P0*S1-P1*S0))
S_23 = eta0123/m0*((H02*H31-H03*H21)*(P2*S3-P3*S2)                    &
+ (H02*H32-H03*H22)*(P3*S1-P1*S3) + (H02*H33-H03*H23)*(P1*S2-P2*S1)   & 
+ (H12*H32-H13*H22)*(P0*S3-P3*S0) + (H12*H33-H13*H23)*(P2*S0-P0*S2)   &
+ (H22*H33-H23*H23)*(P0*S1-P1*S0))

END SUBROUTINE calc_spin_antisym





SUBROUTINE calc_riemann(r,theta)
IMPLICIT NONE
real(kind = dp) :: r,theta,M,Dl,Sg


M = 1.0_dp
Sg = r**2+a**2*cos(theta)**2  
Dl = r**2 + a**2 - 2.0_dp*M*r

R_0101 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*Dl+a**2*sin(theta)**2)/(Sg**3*Dl)


R_0102 = 3.0_dp*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**3


R_0103 = 0.0_dp


R_0112 = 0.0_dp


R_0113 = -M*r*a*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2  &
/(Sg**3*Dl)


R_0123 = M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)           &
*sin(theta)*cos(theta)/Sg**3


R_0202 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(Dl+2.0_dp*a**2*sin(theta)**2)/Sg**3


R_0203 = 0.0_dp


R_0212 = 0.0_dp


R_0213 = M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)             &
*sin(theta)*cos(theta)/Sg**3


R_0223 = M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**3


R_0303 = M*r*Dl*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**3


R_0312 = -M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**2


R_0313 = 0.0_dp


R_0323 = 0.0_dp


R_1212 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg*Dl)


R_1213 = 0.0_dp


R_1223 = 0.0_dp


R_1313 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*((r**2+a**2)**2+2.0_dp*a**2*Dl*sin(theta)**2)  &
/(Sg**3*Dl)


R_1323 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)**3*cos(theta)/Sg**3


R_2323 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*Dl*sin(theta)**2)  &
/Sg**3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

R_0201 = R_0102


R_0301 = R_0103


R_1201 = R_0112


R_1301 = R_0113


R_2301 = R_0123


R_0302 = R_0203


R_1202 = R_0212


R_1302 = R_0213


R_2302 = R_0223


R_1203 = R_0312


R_1303 = R_0313


R_2303 = R_0323


R_1312 = R_1213


R_2312 = R_1223


R_2313 = R_1323


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R0_001 = 0.0_dp


R0_002 = 0.0_dp


R0_003 = 2.0_dp*sin(theta)**2*M**2*a*r**2*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**4


R0_012 = -2.0_dp*(3.0_dp*r**2-a**2*cos(theta)**2)*r*M**2*cos(theta)*a**2*sin(theta)/(Sg**3*Dl)


R0_013 = 0.0_dp


R0_023 = 0.0_dp


R0_101 = M*r*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg**3*Dl)


R0_102 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2*M*r)*M*a**2*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_103 = 0.0_dp


R0_112 = 0.0_dp


R0_113 = 3.0_dp*sin(theta)**2*a*M*r*(r**2+a**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg**3*Dl)


R0_123 = -cos(theta)*sin(theta)*M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*sin(theta)**2*Dl)/(Sg**3*Dl)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R0_201 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4*M*r)*a**2*sin(theta)*cos(theta)*M/(Sg**3*Dl)


R0_202 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(r**2+a**2+2*a**2*sin(theta)**2)/Sg**3


R0_203 = 0.0_dp


R0_212 = 0.0_dp


R0_213 = -(3.0_dp*r**2-a**2*cos(theta)**2)*((a**2+r**2)**2+2.0_dp*a**2*sin(theta)**2*Dl)*M*a*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_223 = -3.0_dp*sin(theta)**2*M*r*a*(a**2+r**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R0_301 = 0.0_dp


R0_302 = 0.0_dp


R0_303 = -sin(theta)**2*M*r*(r**2+3*a**2*sin(theta)**2-3*a**2)*((a**2+r**2)**2-a**2*sin(theta)**2*Dl)/Sg**4


R0_312 = (3.0_dp*r**2-a**2*cos(theta)**2)*((a**2+r**2)**2-a**2*sin(theta)**2*Dl)*M*a*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_313 = 0.0_dp


R0_323 = 0.0_dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R1_001 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(a**2*sin(theta)**2+2.0_dp*Dl)/Sg**4


R1_002 = -3.0_dp*Dl*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**4


R1_003 = 0.0_dp


R1_012 = 0.0_dp


R1_013 = M*r*a*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**4


R1_023 = -a*M*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl*cos(theta)*sin(theta)/Sg**4


R1_101 = 0.0_dp


R1_102 = 0.0_dp


R1_103 = 0.0_dp


R1_112 = 0.0_dp


R1_113 = 0.0_dp


R1_123 = 0.0_dp


R1_201 = 0.0_dp


R1_202 = 0.0_dp


R1_203 = -cos(theta)*sin(theta)*a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl/Sg**3


R1_212 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**2


R1_213 = 0.0_dp


R1_223 = 0.0_dp


R1_301 = -sin(theta)**2*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)/Sg**4


R1_302 = a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)*Dl*cos(theta)*sin(theta)/Sg**4


R1_303 = 0.0_dp


R1_312 = 0.0_dp


R1_313 = -sin(theta)**2*M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*((r**2+a**2)**2+2.0_dp*a**2*Dl*sin(theta)**2)/Sg**4


R1_323 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl*cos(theta)*sin(theta)**3/Sg**4 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R2_001 = -3.0_dp*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**4


R2_002 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*a**2*sin(theta)**2+Dl)/Sg**4


R2_003 = 0.0_dp


R2_012 = 0.0_dp


R2_013 = -a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)*cos(theta)*sin(theta)/Sg**4


R2_023 = -M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**4


R2_101 = 0.0_dp


R2_102 = 0.0_dp


R2_103 = a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/Sg**3


R2_112 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/(Dl*Sg**2)


R2_113 = 0.0_dp


R2_123 = 0.0_dp


R2_201 = 0.0_dp


R2_202 = 0.0_dp


R2_203 = 0.0_dp


R2_212 = 0.0_dp


R2_213 = 0.0_dp


R2_223 = 0.0_dp


R2_301 = a*M*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/Sg**4


R2_302 = M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**4


R2_303 = 0.0_dp


R2_312 = 0.0_dp


R2_313 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)**3/Sg**4


R2_323 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*Dl*sin(theta)**2)*sin(theta)**2/Sg**4 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R3_001 = 0.0_dp


R3_002 = 0.0_dp


R3_003 = 0.0_dp


R3_012 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(a**2*sin(theta)**2-Dl)*M*a*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_013 = 0.0_dp


R3_023 = 0.0_dp


R3_101 = 3.0_dp*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)/(Dl*Sg**3)


R3_102 = -a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*a**2*sin(theta)**2+Dl)*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_103 = 0.0_dp


R3_112 = 0.0_dp


R3_113 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)/(Dl*Sg**3)


R3_123 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*M*a**2*sin(theta)*cos(theta)/(Dl*Sg**3)


R3_201 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(a**2*sin(theta)**2+2.0_dp*Dl)*M*a*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_202 = -3.0_dp*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R3_203 = 0.0_dp


R3_212 = 0.0_dp


R3_213 = -M*a**2*(3*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*sin(theta)*cos(theta)/(Dl*Sg**3)


R3_223 = -M*r*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R3_301 = 0.0_dp


R3_302 = 0.0_dp


R3_303 = -2.0_dp*M**2*a*r**2*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**4


R3_312 = 2.0_dp*M**2*r*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/(Dl*Sg**3)


R3_313 = 0.0_dp


R3_323 = 0.0_dp

END SUBROUTINE calc_riemann






END MODULE tensors
