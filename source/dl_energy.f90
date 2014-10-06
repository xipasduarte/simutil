! dl_energy
!
! Created on: Jul 30, 2014
! Author: Pedro Duarte
! License: MIT License
! Version: 0.0.1
! Description: Extract data from DL_POLY's OUTPUT file.

program dl_energy
	implicit none
	! Variable declaration
	! Allocatable
	! Allocated
	character :: arg*20, label*9, simul*30, labels(3,10)*9, keys*80
	integer :: narg, i, j, k, ts_tot=0, dts=0, ts=1, curr_ts
	
	! Handle command line arguments
	! Check if any arguments are found
	narg = command_argument_count()
	! Loop over the arguments
	if(narg>0)then
		! loop across options
		do i=1,narg
			call get_command_argument(i,arg)
			select case(adjustl(arg))
				case("--version", "-v")
					write(*,*) "v0.0.1"
					stop
				case("--help","-h")
					write(*,*) "Description and tag list with definition."
					stop
				case default
					write(*,*)"Option unknown: ", adjustl(arg)
					stop
			end select
		end do
	end if
	

	open(8, file="OUTPUT", status="old", action="read") ! Open OUTPUT
	
	do
		if(ts_tot==0 .and. dts==0) then
			read(8,"(A30)", advance="no") simul
			if(simul==" selected number of timesteps ") then
				read(8,*) ts_tot
			else if(simul==" data printing interval       ")then
				read(8,*) dts
			else
				read(8,*)
			end if
		else
			read(8,"(A9)") label
			
			if(label == "---------") then
				read(8,"(A9,9(3X,A9))") labels(1,1), labels(1,2), labels(1,3), labels(1,4), labels(1,5),&
				labels(1,6), labels(1,7), labels(1,8), labels(1,9), labels(1,10)
				read(8,"(A9,9(3X,A9))") labels(2,1), labels(2,2), labels(2,3), labels(2,4), labels(2,5),&
				labels(2,6), labels(2,7), labels(2,8), labels(2,9), labels(2,10)
				read(8,"(A9,9(3X,A9))") labels(3,1), labels(3,2), labels(3,3), labels(3,4), labels(3,5),&
				labels(3,6), labels(3,7), labels(3,8), labels(3,9), labels(3,10)
			
				exit
			end if
		end if
	end do
	
	! Ask for the keys
	do i=1,3
		do j=2,10			
			if(j==10) then
				write(*,"(I3,A1,A9)") (i-1)*10+j,":",labels(i,j)
			else
				write(*,"(I3,A1,A9)", advance="no") (i-1)*10+j,":",labels(i,j)
			end if
		end do
	end do
	
	write(*,*) "Keys (space separation): "
	read(*,*) keys
	
	! Discard empty line
	read(8,*)
	
	k = 1
	do while(k<=8)
		read(8,"(I9)") curr_ts
		
		
					
end program dl_energy