! surfangle
!
! Created on: Apr 28, 2014
! Author: Pedro Duarte and Miguel Teixeira
! License: MIT License
! Description: Calculate the angle from the plane on the base of the 
! molecule, to a mean segment representing the molecule. Uses DL_POLY
! HISTORY files as input.

program surfangle

	implicit none
	! Variable declaration
	character(len=80) :: fmt, label
	character (len=132) :: skip
	integer :: index, nmol, status, frame, global_nmol
	real :: x1, y1, z1, x2, y2, z2, global_angle, angles, angle, a, c
	real, parameter :: pi = 4 * atan(1.0)
	logical :: start

	open(8, file="HISTORY", status="old", action="read") ! To do the actual conversion
	open(9, file="HISTORY-angle", status="replace", action="write") ! File to write output
	open(10, file="HISTORY-angles", status="replace", action="write") ! File to write output
	
	!File HEADER
	write(9,*) "Timestep average angles to surface (Degrees). Average over all simulation at the end of file."
	write(10,*) "Angles to surface (Degrees)."

	!File BODY
	start = .TRUE.
	global_angle = 0
	global_nmol = 0
	
	angles = 0
	nmol = 0
	
	! Discard HEADER
	read(8,*) skip; read(8,*) skip;
	
	do
		read(8,*,iostat=status) label, index
		
		if(status .eq. -1) then
			EXIT
		end if
	
		if(label=="timestep") then
			if(start) then
				start = .FALSE.
				
				write(10,*) "timestep    ", index
			else
				write(10,*) "timestep    ", index
				
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
			
			angle = acos(a/c)
			angles = angles + angle
			nmol = nmol + 1
			
			! Write all angles
			write(10,*) index, angle
		else
			read(8,*,iostat=status) skip
		end if
		
		if(status .eq. -1) then
			EXIT
		end if
	end do
	
	! Write final angle
	write(9,*) frame, (angles/nmol*360/2/pi)

	! Write Average Angle
	write(9,*) "Average Angle", global_angle/global_nmol*360/2/pi
	
	close(8)
	close(9)

end program surfangle
