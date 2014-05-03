! surfplane
!
! Created on: Apr 30, 2014
! Author: Pedro Duarte and Miguel Teixeira
! License: MIT License
! Description: Determine the surface cluster of water molecules on the surface
! in DL_POLY simulations.

program surfplane

	implicit none
	! Variable declaration
	integer :: n
	parameter ( n=4000 )
	real, dimension(n*5) :: plane
	integer :: nmol, status, i, j, k, index
	character :: label*8, skip*80
	real :: x, y, z, rW, dist, maxz, viz
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-plane", status="replace", action="write") ! File to write output
	
	nmol = 0
	rW = 1.6
	status = 0
	maxz = 0
	
	! Removes the first two header lines
	read(8,*) skip
	read(8,*) skip
	
	! Count how many neighbors each point has
	do	
		read(8,*,iostat=status) label, index
		
		if(status==-1) then
			EXIT
		end if
		
		if(label=="timestep") then ! Discard timestep and box information
			read(8,*) skip
			read(8,*) skip
			read(8,*) skip
			
		else if(label=="OHW") then ! Handle water molecule
			nmol = nmol + 1
			
			read(8,*) x, y, z
			
			if(z>maxz) then
				maxz = z
			end if
			
			j = (nmol-1)*(5)
			
			do i=1,(nmol-1)
				
				dist = sqrt( (plane((i-1)*5+3)-x)**2 + (plane((i-1)*5+4)-y)**2 + (plane((i-1)*5+5)-z)**2 )
				if(dist < 4*rW) then
					plane((i-1)*5+2) = plane((i-1)*5+2) + 1
				end if
				
			end do
			
			plane(j+1) = index
			plane(j+2) = 0
			plane(j+3) = x
			plane(j+4) = y
			plane(j+5) = z
			
		else
			read(8,*,iostat=status) skip
		end if
		
		if(status==-1) then
			EXIT
		end if
	end do
	
	close(8)
	
	! Calculate average neighbors
	viz = 0
	do i=2,n*5,5
		viz = viz + plane(i)
	end do
	viz = viz/n
	
	! Determine the surface points
	do i=1,n
		
		if(plane((i-1)*5+2) < viz) then
			if( (maxz-plane((i-1)*5+5)) < 2*rW ) then
				write(9,*) plane((i-1)*5+1), plane((i-1)*5+2), plane((i-1)*5+5)
			end if
		end if
		
	end do
	
	close(9)
	
end program surfplane