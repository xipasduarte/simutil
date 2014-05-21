! cutz
!
! Created on: May 21, 2014
! Author: Pedro Duarte and Miguel Teixeira
! License: MIT License
! Description: Remove all molecules below a defined z value.

program cutz

	implicit none
	! Variable declaration
	character :: label*8, fmt1*30, fmt3*30, fmt2*30, fmt4*30
	character (len=20) :: skip
	integer :: index, i, status, a, b, c
	real :: x, y, z, vx, vy, vz, ax, ay, az, minz, xbox, ybox, zbox, d

	open(8, file="CONFIG", status="old", action="read") ! Open CONFIG
	open(9, file="CONFIG-cutz", status="replace", action="write") ! Create output CONFIG
	
	! Call for z value
	write(*,*) "Provide the value of z from which to cut."
	read(*,*) minz
	
	! Variable inicialization
	fmt1 = '(3I10.1, 4X, ES16.10E2)'
	fmt2 = '(3(F20.12))'
	fmt3 = '(3(f20.10))'
	fmt4 = '(A8, I10)'
	
	index = 1 
	
	! File HEADER
	read(8,*) skip
	write(9,*) skip
	
	read(8,*) a, b, c, d
	write(9,fmt1) a, b, c, d
	
	do i=1,3
		read(8,*) xbox, ybox, zbox
		write(9,fmt2) xbox, ybox, zbox
	end do
	
	! Write molecules above z
	do
		read(8,*,iostat=status) label
		if(status .eq. -1) then
			EXIT
		end if
	
		if(label=="HOW") then
			read(8,*) x,y,z
			read(8,*) vx, vy, vz
			read(8,*) ax, ay, az
			
			if(z .gt. minz) then
				write(9,fmt4) label, index
				index = index + 1
				write(9,fmt3) x, y, z
				write(9,fmt3) vx, vy, vz
				write(9,fmt3) ax, ay, az
				
				do i=1,2
					read(8,*) label
					write(9,fmt4) label, index
					index = index + 1
					
					read(8,*) x, y, z
					read(8,*) vx, vy, vz
					read(8,*) ax, ay, az
					write(9,fmt3) x, y, z
					write(9,fmt3) vx, vy, vz
					write(9,fmt3) ax, ay, az
				end do
				
			else
				do i=1,8
					read(8,*) skip
				end do
			end if
		else
			write(9,fmt4) label, index
			index = index + 1
			read(8,*) x, y, z
			read(8,*) vx, vy, vz
			read(8,*) ax, ay, az
			write(9,fmt3) x, y, z
			write(9,fmt3) vx, vy, vz
			write(9,fmt3) ax, ay, az
		end if
	end do
	
	close(8); close(9);

end program cutz
