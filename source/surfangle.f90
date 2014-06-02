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
	character :: fmt*6='(f6.3)', label*8, top*8, bottom*8
	character (len=132) :: skip
	integer :: index, nmol=0, status, frame, global_nmol=0, hist_slot, i
	integer, dimension(9) :: hist
	real :: x1, y1, z1, x2, y2, z2, global_angle=0., angles=0., angle, ori_angle, a, c, d, xbox, ybox
	real, parameter :: pi = 4 * atan(1.0)
	logical :: start=.true.

	open(8, file="HISTORY", status="old", action="read") ! Open HISTORY
	open(9, file="HISTORY-angle", status="replace", action="write") ! Create and open output for average angles
	open(10, file="HISTORY-angles", status="replace", action="write") ! Create and open output for all angles
	open(11, file="HISTORY-hist", status="replace", action="write") ! Create and open output for histograms
	
	! File HEADER
	write(9,*) "Timestep average angles to surface (Degrees). Average over all simulation at the end of file."
	write(10,*) "Angles to surface (Degrees)."
	!write(11,*) "ts          angle       frequency   "
	
	! Variable initialization
	top = "CT"
	bottom = "OH"
	hist=0
	
	! Discard HEADER
	read(8,*) skip; read(8,*) skip;
	
	do
		read(8,*,iostat=status) label, index
		
		if(status .eq. -1) then
			EXIT
		end if
	
		if(label=="timestep") then
			if(start) then
				start = .false.
				
				read(8,*) xbox
				read(8,*) skip, ybox
				read(8,*)
				
				write(10,*) "timestep    ", index
			else
				write(10,*) "timestep    ", index
				
				! Write average angle for previous timestep
				write(9,*) frame, angles/nmol
			
				! Global Average Angle
				global_angle = global_angle + angles
				global_nmol = global_nmol + nmol
			
				! Write histogram
				do i=1,size(hist)
					write(11,*) frame, i*10, hist(i)
				end do
			
				! Reset variables
				angles = 0
				nmol = 0
				hist = 0
			
				read(8,*)
				read(8,*)
				read(8,*)
			end if
			
			frame = index
	
		else if(label == bottom) then
			read(8,'((2X,ES11.4E3),2(1X,ES11.4E3))') x1, y1, z1
	
		else if(label == top) then
			read(8,'((2X,ES11.4E3),2(1X,ES11.4E3))') x2, y2, z2
		
			a = sqrt( (x2-x1)**2 + (y2-y1)**2 )
			
			! Boundary conditions correction
			if(abs(x1-x2) > (xbox/2)) then
				write(*,*) x1, x2, xbox
				if( x2 < 0 ) then
					x2 = x2 + xbox
				else          
					x2 = x2 - xbox
				end if
				
			end if
			
			if(abs(y1-y2) > (ybox/2)) then
				if( y2 < 0 ) then
					y2 = y2 + ybox
				else
					y2 = y2 - ybox
				end if
				
			end if
			
			c = sqrt( (x2-x1)**2 + (y2-y1)**2 +(z2-z1)**2 )
			
			! Calculate orientation angle
			d = abs(x2-x1)
			
			ori_angle = acos(d/c)*360/2/pi
			if (x2 > x1 .and. y2 < y1) then
				ori_angle = 360 - ori_angle
			end if
			if (x2 < x1) then
				if(y2 > y1) then
					ori_angle = 180 - ori_angle
				end if
				if(y2 < y1) then
					ori_angle = 180 + ori_angle
				end if
			end if			
			
			angle = acos(a/c)*360/2/pi
			angles = angles + angle
			nmol = nmol + 1
			
			! Write to all angles
			write(10,*) index, angle, ori_angle
			
			! Create histogram
			hist_slot = int(floor(angle/10))
			hist(hist_slot) = hist(hist_slot) + 1 
		else
			read(8,*,iostat=status)
		end if
		
		if(status .eq. -1) then
			EXIT
		end if
	end do
	
	! Write final angle and histogram
	write(9,*) frame, angles/nmol
	do i=1,size(hist)
		write(11,*) frame, i*10, hist(i)
	end do

	! Write Average Angle
	write(9,*) "Average Angle", global_angle/global_nmol
	
	close(8); close(9); close(10);

end program surfangle
