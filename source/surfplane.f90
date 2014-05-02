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
	integer, dimension(304000) :: plane
	integer :: nmol, status, i, j, k, index
	character :: label
	real :: x, y, z, rW, dist
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-plane", status="replace", action="write") ! File to write output
	
	nmol = 0
	rW = 1.6
	
	do
		read(8,*,iostat=status) label, index
		
		if(status .ne. 0 .or. status .ne. 5010) then
			EXIT
		end if
		
		if(label=="OHW") then
			nmol = nmol + 1
			
			read(8,*) x, y, z
			
			j = (nmol-1)*(3*25+1)
			
			do i=1,j,(3*25+1)
				
				dist = sqrt( (plane(i+1)-x)**2 + (plane(i+2)-y)**2) + (plane(i+3)-z)**2 )
				if(dist < 4*rW) then
					do k=i+4,i*(3*25+1),3
						if(ABS(plane(k)) + ABS(plane(k+1)) + ABS(plane(k+2)) == 0) then
							plane(k) = x
							plane(k+1) = y
							plane(k+2) = z
						end if
					end do
				end if
				
			end do
			
			plane(j+1) = index
			plane(j+2) = x
			plane(j+3) = y
			plane(j+3) = z
		end if
	end do
	
	do i=1,304000,(3*25+1)
		
		index = plane(i)
		
		if(index==0) then
			EXIT
		end if
		
		if(ABS(plane(i*(3*25+1)-2)) + ABS(plane(i*(3*25+1)-1)) + ABS(plane(i*(3*25+1))) == 0) then
			write(9,*) index
		end if
		
	end do
	
end program surfplane
