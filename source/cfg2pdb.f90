! cfg2pdb
! 
! Revised on: May 1, 2014
! Author: Zé Nuno
! Revision Author: Pedro Duarte
! License: MIT License
! Description: Converts a DL_POLY configuration file 
! into a rasmol readable pdb file.

program cfg2pdb                                           

	implicit none

	character :: header*80, atmnam*4, fmt_header*20
	real :: x, y, z
	integer :: nmol, levcnf, imcon, status, i

	open(2, file='CONFIG', status="old", action="read")
	open(3, file='RASMOL', status="replace", action="write")

	! Variable initialization
	fmt_header = '(A80,/,2(I10))'

	read(2,fmt_header) header, levcnf, imcon
	write(3,*) header
	
	read(2,'(2/)')
	
	do
		read (2,*,iostat=status) atmnam, nmol
		
		if(status .ne. 0) then
			EXIT
		end if
		read(2,'(3F20.0)') x, y, z

		if(len_trim(atmnam)<4) then
			write(3,'(A6,I5,2X,A3,14X,3F8.3,23X,A1)') &
			'ATOM  ', nmol, atmnam, x, y, z, atmnam
		else
			write(3,'(A6,I5,1X,A4,14X,3F8.3,23X,A1)') &
			'ATOM  ', nmol, atmnam, x, y, z, atmnam
		end if
		
        do i=1,levcnf
			read(2,*)
		end do
	end do

	close (2) ; close (3)

end program

