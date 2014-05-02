! cfg2pdb
! 
! Revised on: May 1, 2014
! Author: Zé Nuno
! License: MIT License
! Description: Converts a DL_POLY configuration file 
! into a rasmol readable pdb file.

program cfg2pdb                                           

	implicit none

	character header*80
	character atmnam*8, end*3
	real(8) xxx, yyy, zzz
	integer nmol, levcnf, imcon

	open (2, file = 'CONFIG') ; open (3, file = 'RASMOL')



	read (2,'(A80,/,2(i10))') header, levcnf, imcon
	if (imcon.gt.0) read (2,'(2(/))')
	write (3,*) header ; write (3,*) levcnf, imcon



	nmol = 1

        read (2,'(A8)') atmnam

	do while (atmnam.ne."END     ")

		read (2,'(3(f20.0))') xxx, yyy, zzz

		if (levcnf.gt.0) read (2,'(/)')

		write (3,'(A6,I5,2X,A3,14X,3F8.3, I3)') &
		'ATOM  ', nmol, atmnam, xxx, yyy, zzz
		nmol = nmol + 1

                read (2,'(A8)') atmnam

	end do

	close (2) ; close (3)

end program

