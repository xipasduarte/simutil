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
	parameter ( n=12 )
	real, dimension(4000*(3*n+1)) :: plane
	integer :: nmol, status, i, j, k, index
	character :: label*8, skip*80
	real :: x, y, z, rW, dist, maxz
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-plane", status="replace", action="write") ! File to write output
	
	nmol = 0
	rW = 1.6
	status = 0
	maxz = 0
	
	! Removes the first two header lines
	read(8,*) skip
	read(8,*) skip
	
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
			
			j = (nmol-1)*(3*n+1)
			
			do i=1,j,(3*n+1)
				
				dist = sqrt( (plane(i+1)-x)**2 + (plane(i+2)-y)**2 + (plane(i+3)-z)**2 )
				if(dist < 4*rW) then
					do k=i+4,i+(3*n+1),3
						if(ABS(plane(k)) + ABS(plane(k+1)) + ABS(plane(k+2)) == 0) then
							plane(k) = x
							plane(k+1) = y
							plane(k+2) = z
							
							EXIT
						end if
					end do
				end if
				
			end do
			
			plane(j+1) = index
			plane(j+2) = x
			plane(j+3) = y
			plane(j+4) = z
			
		else
			read(8,*,iostat=status) skip
		end if
		
		if(status==-1) then
			EXIT
		end if
	end do
	
	close(8)
	
	do i=1,4000
		
		if(ABS(plane(i*(3*n+1)-2)) + ABS(plane(i*(3*n+1)-1)) + ABS(plane(i*(3*n+1))) == 0) then
			if( (maxz-plane((i-1)*(3*n+1)+4)) < 2*rW ) then
				write(9,*) plane((i-1)*(3*n+1)+1), plane((i-1)*(3*n+1)+4)
			end if
		end if
		
	end do
	
	close(9)
	
end program surfplane