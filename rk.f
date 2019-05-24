MODULE runge_kutta
USE SHARED_CONSTANTS
USE Dixon_equations
USE spin_condition 
USE tensors
IMPLICIT NONE

CONTAINS

SUBROUTINE rk(ystart,yy,nstep,t_final)
IMPLICIT NONE 
real(kind=dp), DIMENSION(:) :: ystart
real(kind=dp), DIMENSION(3) :: COM, COR !Centre of mass, Centre of radiation
real(kind=dp), DIMENSION(size(ystart)) :: y,y1, dy1, y2, y3, dy2, dy3,&
dy4,dy5,dy6  
REAL(KIND=dp), DIMENSION(:,:), INTENT(OUT) :: yy
INTEGER(kind=dp) :: k,nstep, finalstep,trigger, step_number, j
real(kind=dp) :: t_final
real(kind=dp) :: mm, xOUT,yOUT,zOUT
real(kind=dp) :: start_time, end_time
real(kind=dp) :: Sx, Sy, Sz, thetaSL, phiSL, psi
real(kind=dp) :: Bx,By,Bz,Cx,Cy,Cz, blah
real(kind=dp), dimension(:,:), allocatable :: dataset, AllData, Cartesian
real(kind=dp), dimension(:,:), allocatable :: INERTIA, ST1, ST2, ST3,MT1, MT2,MT3
real(kind=dp) :: nx, ny, nz, nr, ntheta,nphi,t0
real(kind=dp) :: Ry_x, Ry_y, Ry_z, Rz_x, Rz_y,Rz_z
real(kind=dp) :: psi0, phi_start
real(kind=dp) :: start_phi,magnitude
real(kind=dp) :: zderiv
real(kind=dp) :: tau, signR, signR_new
real(kind=dp) :: phi_old, phi_new
real(kind=dp) :: periapsis_interpolated, periapsis_averaged
real(kind=dp) :: FinalPhi, dphi, TargetPhi
integer(kind=dp) :: TargetJ,MetaK,PlotK,PlotRows, InertRows, N, Nc
real(kind=dp) :: hplus, hcross
real(kind=dp),dimension(6) :: harray
real(kind=dp), DIMENSION(size(ystart)):: k1,k2,k3,k4,k5,k6,yscal,ratio
real(kind=dp), DIMENSION(size(ystart)):: y4,y5,y6,ynew,yerr,DeltaErr
real(kind=dp) :: errmax
real(kind=dp), dimension(:,:), allocatable :: ITensor, hOUT
real(kind=dp), dimension(:,:), allocatable :: IDerivs, SDerivs1, Sderivs2, SDerivs3
real(kind=dp), dimension(:,:), allocatable :: MDerivs1, Mderivs2, MDerivs3
real(kind=dp), dimension(:,:), allocatable :: MTot, STot,TTot
real(kind=dp), dimension(:), allocatable :: PlayIn, PlayOut
integer(kind=dp) :: order
!Create massive array for storing output
PlotRows = 1e8
ALLOCATE(AllData(PlotRows, 8)) !4x BL, 4* 4 vel









!Initialise counter
PlotK = 1



!Create plot save file
open(14,file=PLOToutfile,status = 'replace')
close(14)






! Setup and Initialize
y(:)=ystart(:)                              
k = 0
MetaK = 0
!CALL derivs(y, dy1)

!Define phi variables
FinalPhi = 2.0_dp*PI*N_orbit
phi_start =  ystart(4)

!Printeger = 1

11 do while ( abs(y(4) - phi_start) .LT. FinalPhi) 
 
   !Calculate k1 
   CALL derivs(y, dy1)
   k1(:) = h*dy1(:)
   
   
   
   !Caclulate k2
   y2(:) = y + B21*k1
   CALL derivs(y2, dy2)
   k2(:) = h*dy2(:)
   

   !Caclulate k3
   y3(:) = y + B31*k1 + B32*k2
   CALL derivs(y3, dy3)
   k3(:) = h*dy3(:)



   !Caclulate k4
   y4(:) = y + B41*k1 + B42*k2 + B43*k3
   CALL derivs(y4, dy4)
   k4(:) = h*dy4(:)



   !Caclulate k5
   y5(:) = y + B51*k1 + B52*k2 + B53*k3 + B54*k4
   CALL derivs(y5, dy5)
   k5(:) = h*dy5(:)


   !Caclulate k6
   y6(:) = y + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
   CALL derivs(y6, dy6)
   k6(:) = h*dy6(:)

 
   ynew = y + c1*k1  + c3*k3 + c4*k4  +c6*k6
   yerr = y + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6

   DeltaErr = abs(ynew - yerr)
   yscal = abs(y) + abs(k1) + 1.0D-3

   ratio = DeltaERR/yscal
   errmax = escal * maxval(ratio)


   if (errmax .GT. 1.0_dp) then

   !This is not good. Adjust the stepsize and try again
   call ShrinkStepsize(errmax)
   goto 11
   else
   !This is good. Try to increase the stepsize a little bit
   call GrowStepsize(errmax)
   endif



    !Check it doesnt fall below event horizon
    if (y(2) .LT. Rhor) then
    print *, 'Whoops! Fell below event horizon'
    print *, 'STOP: ', y(2), Rhor
    EXIT
    endif


   !Update for the next timestep
   y = ynew

CALL derivs(y, dy1)
if (PlotK .EQ. 1) then

 !Calculate hplus, hcross time series
 call GravWaves(y(2),y(3),y(4), &
                   dy1(2), dy1(3),dy1(4), &
                   dy1(6)/m0, dy1(7)/m0,dy1(8)/m0, &
                   hplus,hcross,harray)
   
endif

   AllData(PlotK,1:4) = y(1:4) !BL coordinates
   !AllData(PlotK,5:8) = y(5:8)/m0 !$ momentum/m = 4 velocity
   AllData(PlotK,5:8) = dy1(1:4) !4 velocity


!    !For plotting the trajectory
     mm = sqrt(y(2)**2 + a**2)
     xOUT = mm*sin(y(3))*cos(y(4))
     yOUT = mm*sin(y(3))*sin(y(4))
     zOUT = mm*cos(y(3))




   !  PlotDataSet(PlotK,1) = xOUT
   !  PlotDataSet(PlotK,2) = !yOUT
  !   PlotDataSet(PlotK,3) = zOUT
     
     !PlotDataSet(PlotK,4) = harray(1) !* OBSR/m0
     !PlotDataSet(PlotK,5) = harray(2)

     if (PlotK .EQ. 1) then
     !print *, y(1),m0*xOUT*xOUT, m0*yOUT*yOUT,m0*zOUT*zOUT, m0*xOUT*yOUT, m0*xOUT*zOUT, m0*yOUT*zOUT
     continue
     endif
     !INERTIA(PlotK,1) = y(1)          !t
     !INERTIA(PlotK,2) = m0*xOUT*xOUT  !Ixx
     !INERTIA(PlotK,3) = m0*yOUT*yOUT  !Iyy
     !INERTIA(PlotK,4) = m0*zOUT*zOUT  !Izz
     !INERTIA(PlotK,5) = m0*xOUT*yOUT  !Ixy
     !INERTIA(PlotK,6) = m0*xOUT*zOUT  !Ixz
     !INERTIA(PlotK,7) = m0*yOUT*zOUT  !Iyz
!


     PlotK = PlotK  + 1
     if (PlotK .GT. PlotRows) then
     print *, 'Need a bigger array to save output'
     endif



02 enddo











Nvector(1) = a*sin(OBSTheta)*cos(OBSPhi)
Nvector(2) = a*sin(OBSTheta)*sin(OBSPhi)
Nvector(3) = cos(OBSTheta)



!print *, AllData(1,:)


!Set Number of entries
N = PlotK - 1


!Calculate cartesian components
ALLOCATE(Cartesian(N,6))

Cartesian(:,1) = sqrt(AllData(1:N,2)**2 + a**2)*sin(AllData(1:N,3))*cos(AllData(1:N,4)) !x
Cartesian(:,2) = sqrt(AllData(1:N,2)**2 + a**2)*sin(AllData(1:N,3))*sin(AllData(1:N,4)) !y
Cartesian(:,3) = sqrt(AllData(1:N,2)**2 + a**2)*cos(AllData(1:N,3)) !z

Cartesian(:,4) = AllData(1:N,2)*AllData(1:N,6)*sin(AllData(1:N,3))*cos(AllData(1:N,4)) / (sqrt(AllData(1:N,2)**2 + a**2)) + &
                 sqrt(AllData(1:N,2)**2 + a**2) * cos(AllData(1:N,4))*cos(AllData(1:N,3))*AllData(1:N,7) -&
                 sqrt(AllData(1:N,2)**2 + a**2) * sin(AllData(1:N,3))*sin(AllData(1:N,4))*AllData(1:N,8) !vx
                 

Cartesian(:,5) = AllData(1:N,2)*AllData(1:N,6)*sin(AllData(1:N,3))*sin(AllData(1:N,4)) / (sqrt(AllData(1:N,2)**2 + a**2)) + &
                 sqrt(AllData(1:N,2)**2 + a**2) * sin(AllData(1:N,4))*cos(AllData(1:N,3))*AllData(1:N,7) +&
                 sqrt(AllData(1:N,2)**2 + a**2) * sin(AllData(1:N,3))*sin(AllData(1:N,4))*AllData(1:N,8) !vy



Cartesian(:,6) = AllData(1:N,6) * cos(AllData(1:N,3)) - AllData(1:N,2)*sin(AllData(1:N,3))*AllData(1:N,7) !vz











!Calculate the mass quadrupole I, mass octupole M and current Quadrupole S

!!!!!!-------------------------!!!!!!!!!!!! Mass Quadrupole
Nc = 7
ALLOCATE(INERTIA(N,Nc))
INERTIA(:,1) = AllData(1:N,1)
INERTIA(:,2) = m0*Cartesian(:,1)*Cartesian(:,1) !xx
INERTIA(:,3) = m0*Cartesian(:,2)*Cartesian(:,2) !yy
INERTIA(:,4) = m0*Cartesian(:,3)*Cartesian(:,3) !zz
INERTIA(:,5) = m0*Cartesian(:,1)*Cartesian(:,2) !xy
INERTIA(:,6) = m0*Cartesian(:,1)*Cartesian(:,3) !xz
INERTIA(:,7) = m0*Cartesian(:,2)*Cartesian(:,3) !yz

ALLOCATE(IDerivs(N,Nc-1))

order = 2
call FiniteDifference(N,Nc, INERTIA,IDerivs,order)




!!!!!!-------------------------!!!!!!!!!!!! Current Quadrupole

!3 arrays to take product with teh obsvector
ALLOCATE(ST1(N,Nc))
ALLOCATE(ST2(N,Nc))
ALLOCATE(ST3(N,Nc))

ST1 = INERTIA
ST2 = INERTIA
ST3 = INERTIA

!Define the Current quadrupole tensor
do j = 2,Nc
ST1(:,j) = Cartesian(:,4) * ST1(:,j)
ST2(:,j) = Cartesian(:,5) * ST1(:,j)
ST3(:,j) = Cartesian(:,6) * ST1(:,j)
enddo





!Take the secondderivative of the current quadrupole tensor

ALLOCATE(SDerivs1(N,Nc-1))
ALLOCATE(SDerivs2(N,Nc-1))
ALLOCATE(SDerivs3(N,Nc-1))



call FiniteDifference(N,Nc, ST1,SDerivs1,order)
call FiniteDifference(N,Nc, ST2,SDerivs2,order)
call FiniteDifference(N,Nc, ST3,SDerivs3,order)


!Takethe product with the obs vector
ALLOCATE(STot(N,Nc-1))
STot = Nvector(1)*SDerivs1 + Nvector(2)*SDerivs2 + Nvector(3)*SDerivs3




!!!!!!-------------------------!!!!!!!!!!!! Mass Octupole

ALLOCATE(MT1(N,Nc))
ALLOCATE(MT2(N,Nc))
ALLOCATE(MT3(N,Nc))

MT1 = INERTIA
MT2 = INERTIA
MT3 = INERTIA


!Define the mass octupole tensor
do j = 2,Nc
MT1(:,j) = Cartesian(:,1) * MT1(:,j)
MT2(:,j) = Cartesian(:,2) * MT1(:,j)
MT3(:,j) = Cartesian(:,3) * MT1(:,j)
enddo




!Take the THIRD derivative
ALLOCATE(MDerivs1(N,Nc-1))
ALLOCATE(MDerivs2(N,Nc-1))
ALLOCATE(MDerivs3(N,Nc-1))

order = 3
call FiniteDifference(N,Nc, MT1,MDerivs1,order)
call FiniteDifference(N,Nc, MT2,MDerivs2,order)
call FiniteDifference(N,Nc, MT3,MDerivs3,order)


!Take the product with the obs vector

ALLOCATE(MTot(N,Nc-1))
MTot = Nvector(1)*MDerivs1 + Nvector(2)*MDerivs2 + Nvector(3)*MDerivs3


!Now calculate bit in the brackets

ALLOCATE(TTot(N,Nc-1))
TTot = IDerivs - 2.0_dp*STot + MTot



!And calculate hplus/hcross strains

ALLOCATE(hOUT(N,2))
call GW2(N,TTot,hOUT)
print *, 'TTOT:', IDerivs(2,1), STot(2,1), MTot(2,1)

print *, 'hout:', hOUT(2,1)*OBSR/m0, hOUT(2,2)*OBSR/m0



!Save output into one file 


!For plotting the trajectory
open (14, file = PLOToutfile, status='replace',action='write')
do j = 1,N
write(14,*) Cartesian(j,1), Cartesian(j,2), Cartesian(j,3), &
            AllData(j,1)/convert_s, hOUT(j,1), hOUT(j,2), &
            hOUT(j,1)*OBSR/m0, hOUT(j,2)*OBSR/m0
enddo
close(14)


print *, 'NUM ROWS  =', N
open(14,file='TEST_INERTIA.txt',status = 'replace', action='write')
do j = 1,N
write(14,*) INERTIA(j,1), hOUT(j,1)*OBSR/m0, hOUT(j,2)*OBSR/m0
enddo
close(14)





!N = 1d4
!print *, INERTIA(1,:)


!Allocate the data arrays
!ALLOCATE(ITensor(N,6))
!call SecondDerivatives(N, INERTIA,ITensor)


!call GW2(N,ITensor,hOUT)


!Write datasets to disk

!PlotDataSet(1,11) = FundFreq
!PlotDataSet(1,13) = timescale


!print *, PlotK

!First write the file for plotting trajectory 
!open(14,file=PLOToutfile,status = 'replace', action='write')
!do j = 1,PlotK-1
!write(14,*) PlotDataSet(j,1), PlotDataSet(j,2),PlotDataSet(j,3),PlotDataSet(j,4), PlotDataSet(j,5), &
 !           PlotDataSet(j,6), PlotDataSet(j,7),PlotDataSet(j,8),PlotDataSet(j,9), PlotDataSet(j,10),&
  !          PlotDataSet(j,11), PlotDataSet(j,12), PlotDataSet(j,13)
!enddo
!close(14)



!open(14,file='TEST_INERTIA.txt',status = 'replace', action='write')
!do j = 1,PlotK-1
!write(14,*) INERTIA(j,1), hOUT(j,1)*OBSR/m0, hOUT(j,2)*OBSR/m0
!enddo
!close(14)


!print *, INERTIA(1,1),INERTIA(1,2)
!open (15, file='TEST_INERTIA.dat', form='unformatted', status='replace')
!write(15) INERTIA
!close(15)


!open (15, file=DATAoutfile, form='unformatted', status='replace')
!write(15) PlotDataSet
!close(15)

!print *, 'FINISH'
!print *, 'Outfiles:'
!print *, PLOToutfile


END SUBROUTINE rk

END MODULE runge_kutta
