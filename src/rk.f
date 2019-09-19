module rungekutta

use parameters
use constants
use derivatives
use GravWaves

implicit none

private


public rk

contains


subroutine rk(y0)

 
!Arguments
real(kind=dp),intent(IN),dimension(entries) :: y0 !initial conditions

!Other
real(kind=dp), dimension(size(y0)) :: y, y1, dy
!an array to store all the data. !Should be faster than dynamically allocating
!However this does require more exploration
!see https://stackoverflow.com/questions/8384406/how-to-increase-array-size-on-the-fly-in-fortran
real(kind=dp), dimension(nrows,size(y0)+3) :: AllData !Big array to save all data
real(kind=dp), dimension(:,:),allocatable :: output, waveforms !smaller array which will be outout
integer(kind=dp) :: i,j !index for saving to array
real(kind=dp) :: mm, xC, yC, zC !Cartesian components

y  = y0


!Save the first row to array
i = 1
AllData(i,1:entries) = y

  

!Integrate
!do while (y(1) .LT. 1.0_dp)
!do while ( abs( y(4) - y0(4)) .LT. N_orbit*2.0_dp*PI)    

!do while (i .LT. 5)

do while (i .LT. 10)

   call RKF(y,y1)
    y = y1
 

    print *, i, y(1), h
   

   !Also calculate the rdot, thetadot, phidot
   call derivs(y,dy)

    !Save the output
    i = i + 1
    AllData(i,1:entries) = y
    AllData(i,entries+1: entries+3) = dy(2:4)


enddo





print *, 'Runge Kutta completed. Start data I/O'
!!!!!!!!!! - Save the output for analysis and plotting - !!!!!!!


!Binary format. See discussion at https://stackoverflow.com/questions/24395686/best-way-to-write-a-large-array-to-file-in-fortran-text-vs-other



!First reallocate to create a smaller array
allocate(output(i,entries+3))
allocate(waveforms(i,2))

output = AllData(1:i, :)


    call calc_GW(i,output,waveforms)






!Now save
open(unit=10,file=BinaryData,status='replace',form='unformatted')
write(10) output
close(10)




if (plot .EQ. 1) then
!Save formatted data for plotting
    open(unit=20,file=PlotData,status='replace',form='formatted')
    do j = 1,i
    if (mod(real(j), coarse) .EQ. 0.0_dp) then
    mm = sqrt(output(j,2)**2 + a**2)
    xC = mm * sin(output(j,3)) * cos(output(j,4))
    yC = mm * sin(output(j,3)) * sin(output(j,4))
    zC = mm * cos(output(j,3)) 
    write(20, *) output(j,1), xC, yC, zC, waveforms(j,1), waveforms(j,2)
    endif
    enddo
    close(20)



endif



!Closing messages
print *, 'Data saved to ', trim(adjustl(BinaryData))
print *, 'Estimated file size is: ', real(i)*real(entries)*real(dp)/1.0d6, ' MB'


end subroutine rk


end module rungekutta
