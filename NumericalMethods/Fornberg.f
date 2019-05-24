Program fornberg

  implicit none
integer, parameter :: dp = selected_real_kind(32,307)
integer, parameter :: Nrows =1d4 
integer (kind=dp) :: j, counts,k, N, i, mn, l
real(kind=dp), dimension(Nrows,2) :: SomeData
real(kind=dp), dimension(:,:), allocatable :: INERTIA, carray
real(kind=dp), dimension(:), allocatable :: xdata,ydata,xtgt
real(kind=dp) :: c1,c2,c3,c4,c5, z
integer(kind=dp),parameter :: nd = 5 !Number of grid points
integer(kind=dp),parameter :: nn = nd - 1 
integer(kind=dp),parameter :: m = 2
real(kind=dp) :: Iout
!Load the data

!Load all of the data
open(unit = 10, file='../TEST_INERTIA.dat', form = 'unformatted')
read(10) SomeData
close(10)


!Get number of rows in dataset
do j=1,Nrows
if (SomeData(j,1) .EQ. 0.0_dp) then
N = j-1
EXIT
endif
enddo



!Allocate the data arrays

ALLOCATE(xdata(0:nd-1))
ALLOCATE(ydata(0:nd-1))
do j = 1,N
!For every point, calculate the second derivative
z =SomeData(j,1)


!print *, j,'z = ',z



if (j .LT. 1 + nn/2) then

!Forward modelling

do l=1,nd
xdata(l-1) = SomeData(j-1+l,1)
ydata(l-1) = SomeData(j-1+l,2)
enddo



elseif (j + nn/2 .GT. N) then
!Backward modelling


do l=1,nd
xdata(l-1) = SomeData(j-(nd-l),1)
ydata(l-1) = SomeData(j-(nd-l),2)
enddo



else
!Central modelling

do l=1,nd
xdata(l-1) = SomeData(j-1-(nd/2)+l,1)
ydata(l-1) = SomeData(j-1-(nd/2)+l,2)
enddo



endif




!print *, 'DATA IN:'
!print *,'X:', xdata
!print *,'Y:', ydata

call SecondDerivative(nd,m,z,xdata,ydata,Iout)



print *, j, Iout


enddo



!Set up zeroes array





!print *, SUM(carray(:,2)*ydata)

end program fornberg



subroutine SecondDerivative(nd,m,z,x,y,Iout)

implicit none
integer, parameter :: dp = selected_real_kind(32,307)
integer(kind=dp) :: nd, m,i,j,k
integer(kind=dp) :: nn,mn
real(kind=dp) :: z
real(kind=dp), dimension(0:nd-1) :: x,y
real(kind=dp), dimension(0:nd-1,0:m) :: c
real(kind=dp) :: c1,c2,c3,c4,c5,Iout

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

Iout =  SUM(c(:,2)*y)


end subroutine SecondDerivative




