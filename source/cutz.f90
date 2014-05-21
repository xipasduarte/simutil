! cutz
!
! Created on: May 21, 2014
! Author: Pedro Duarte and Miguel Teixeira
! License: MIT License
! Description: Remove all molecules below a defined z value.

program cutz

	implicit none
	! Variable declaration
	character :: label*8, fmt_all*5, fmt1*30, fmt2*30, all*60
	character (len=16) :: x, y, z
	integer :: index, i, j, status
	real :: minz, numbz

	open(8, file="CONFIG", status="old", action="read") ! Open CONFIG
	open(9, file="CONFIG-cutz", status="replace", action="write") ! Create output CONFIG
	
	! Call for z value
	write(*,*) "Provide the value of z from which to cut."
	read(*,*) minz
	
	! Variable inicialization
	fmt_all = '(A60)'
	fmt1 = '(A16, 2(A20))'
	fmt2 = '(A8, I10)'
	
	index = 1 
	
	! File HEADER
	do i=1,5
		read(8,fmt_all) all
		write(9,fmt_all) all
	end do
	
	! Write molecules above z
	do
		read(8,*,iostat=status) label
		if(status .eq. -1) then
			EXIT
		end if
	
		if(label=="HOW") then
			read(8,fmt1) x,y,z
			read(z,*) numbz
			
			if(numbz .gt. minz) then
				write(9,fmt2) label, index
				index = index + 1
				write(9,fmt1) x, y, z
				
				do i=1,2
					read(8,fmt_all) all
					write(9,fmt_all) all
				end do
				
				do i=1,2
					read(8,*) label
					write(9,fmt2) label, index
					index = index + 1
					
					do j=1,3
						read(8,fmt_all) all
						write(9,fmt_all) all
					end do
				end do
				
			else
				do i=1,10
					read(8,*)
				end do
			end if
		else
			write(9,fmt2) label, index
			index = index + 1
			do i=1,3
				read(8,fmt_all) all
				write(9,fmt_all) all
			end do
		end if
	end do
	
	close(8); close(9);

end program cutz
