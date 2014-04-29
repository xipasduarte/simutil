! surfangle
!
! Created on: Apr 28, 2014
! Author: Pedro Duarte e Miguel Teixeira
! License: MIT
! Description: Calculate the angle from the plane on the base of the 
! molecule, to a mean segment representing the molecule.

program surfangle

	implicit none
	! Variable declaration
	character(len=80) :: fmt, label
	character (len=132) :: skip
	integer :: index, nmol, status, frame, global_nmol
	real :: x1, y1, z1, x2, y2, z2, global_angle, angles, a, c
	real, parameter :: pi = 4 * atan(1.0)
	logical :: start

	open(8, file="HISTORY", status="old", action="read") ! To do the actual conversion
	open(9, file="HISTORY-ANGLE", status="replace", action="write") ! File to write output
	
	! Read positions and make CONFIG file
	fmt="(3(f20.12))"
	!File HEADER
	write(9,*) "Angles to surface (Degrees). Average over all simulation at the end of file."

	!File BODY
	start = .TRUE.
	global_angle = 0
	global_nmol = 0
	
	angles = 0
	nmol = 0
	status = 0
	
	do while(status .eq. 0)
		read(8,*,iostat=status) label, index
		
		if(status .ne. 0) then
			! Write Average Angle
			write(9,*) "Average Angle", global_angle/global_nmol*360/2/pi
			EXIT
		end if
		
		if(label=="timestep") then
			if(start) then
				start = .FALSE.
			else
				write(9,*) frame, (angles/nmol*360/2/pi)
				
				! Global Average Angle
				global_angle = global_angle + angles
				global_nmol = global_nmol + nmol
				
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
			read(8,*) skip
		end if
	end do

	close(8)
	close(9)

end program surfangle
