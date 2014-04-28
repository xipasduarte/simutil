! grotoconfig
!
! Created on: Mar 9, 2012
! Author: Pedro Duarte e Pedro Morgado

program surfangle

	implicit none
	! Variable declaration
	character(len=80) :: fmt
	character (len=132) :: skip
	integer :: index, nmol, status
	real :: x1, y1, z1, x2, y2, z2, global_angle, angles, a, c
	logical :: start

	open(8, file="HISTORY", status="old", action="read") ! To do the actual conversion
	open(9, file="HISTORY-ANGLE", status="replace", action="write") ! File to write output
	
	! Read positions and make CONFIG file
	fmt="(3(f20.12))"
	!File HEADER
	write(9,*) "HEADER"

	!File BODY
	start = .TRUE.
	global_angle = 0
	angles = 0
	nmol = 0
	
	do while(status .eq. 0)
		read(8,*) label, index
		
		if(label=="timestep") then
			if(start) then
				start = .FALSE.
			else
				write(9,*) frame, angles/nmol
				
				! Global Average Angle
				global_angle = global_angle
				
				! Reset variables
				angles = 0
				nmol = 0
			end if
			
			frame = index
			read(8,*) skip
			read(8,*) skip
			read(8,*) skip
		
		else if(label == "OH") then
			read(8,*) x1, y1, z1
		
		else if(label == "CT") then
			read(8,*) x2, y2, z2
			
			a = sqrt( (x2-x1)**2 + (y2-y1)**2 )
			c = sqrt( (x2-x1)**2 + (y2-y1)**2 +(z2-z1)**2 )
			
			angles = angles + acos(a/c)
			nmol = nmol + 1
		
		else
			read(8,*, iostat=status) skip
		end if
	end do

	close(8)
	close(9)

end program surfangle
