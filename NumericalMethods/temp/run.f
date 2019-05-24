program run

use fornberg

  implicit none

  integer, parameter :: dp=kind(0.d0)  ! double precision


  real(dp) :: z

  integer, parameter   :: nd = 271 - 1
  integer,parameter    :: m = 2
  real(dp)    :: x(0:nd)
  real(dp)    :: yy(0:nd)
  real(dp)   :: c(0:nd, 0:m)
  real(dp) :: t,I
  integer :: ios,counts


open(unit = 10, file='../../TEST_INERTIA.txt', form = 'formatted',action='read')

ios= 0
counts = 0
do while (ios .EQ. 0)
read(10,*,iostat=ios) t, I

if (ios .EQ. 0) then
print *, counts, t,I, ios
x(counts) = t
counts = counts+1
endif
enddo


close(10)



  !z = 0.0_dp
  !x(0) = -2.0_dp
  !x(1) = -1.0_dp
  !x(2) = 0.0_dp
  !x(3) = +1.0_dp
  !x(4) = +2.0_dp

  z = x(0)

  print *, 'Z = ', z
  print *, 'X = ', x
  print *, 'ND = ', nd 
  print *, 'm = ', m
  print *, 'c = ', c 
  call populate_weights(z,x,nd,m,c)

  print *, c(:,2)

end program run
