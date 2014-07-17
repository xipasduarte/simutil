! hbd
!
! Created on: Jun 2, 2014
! Author: Pedro Duarte and Pedro Morgado
! License: MIT License
! Version: 0.1.1
! Description: Determine the distribution of hydrogen bonds in mixtures – both O and H participate
! in hydrogen bonds; only H; only O; and none (the molecule doesn't make any hydrogen bonds).

program hbd

	implicit none
	! Variable declaration
	integer :: narg, i, j, k, status=0, index, natoms, ts=0
	character :: arg*20, label*8, skip
	real :: box(3)
	! Allocatable
	character :: labels(:,:)*8
	integer :: hbdist(:,:)
	
	! Handle command line arguments
	! Check if any arguments are found
	narg = command_argument_count()
	! Loop over the arguments
	if(narg>0)then
		! loop across options
		do i=1,narg
			call get_command_argument(i,arg)
			select case(adjustl(arg))
				case("--labels")
					if(mod(narg-1,2).ne.0) then
						write(*,*) "--labels must have even number of labels"
						stop
					end if
					allocate(labels((narg-i)/2,2))
					do k=i+1,narg,2
						do j=1,2
							call get_command_argument(k+(j-1), arg)
							labels(k/2,j) = trim(arg)
						end do
					end do
					EXIT
				case("--version", "-v")
					write(*,*) "v0.1.1"
					stop
				case("--help","-h")
					write(*,'(A33,(/))') "hbd (Hydrogen Bond Distribution)"
					write(*,*) "Determine the types of hydrogen bonding in mixtures."
					stop
				case default
					write(*,*)"Option unknown: ", adjustl(arg)
					stop
			end select
		end do
	else
		allocate(labels(1,2))
		labels(1,1:2) = (/ "HOT", "OHT" /)
	end if
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-hbd", status="replace", action="write") ! File to write output
	
	! Removes the first two header lines
	read(8,*)
	read(8,*) skip, skip, natoms
	
	! Initialize hbdist array
	allocate(hbdist(size(labels,1)+1,5))
	hbdist(:,:) = 0
	
	do	
		read(8,*,iostat=status) label, index
		
		if(status==-1) then
			EXIT
		end if
		
		if(label=="timestep") then
			ts = ts + 1
			
			read(8,*) box(1)
			read(8,*) skip, box(2)
			read(8,*) skip, skip, box(3)
			
			call ts_hbd(natoms, labels, size(labels,1), hbdist, box)
		end if
	end do
	
	! Average frequencies
	hb_count = hb_count/ts
	
	write(9,"(5(A12))") "HBonds", "None", "Hydrogen", "Oxygen", "Both" 
	write(9,"(A12,4(f12.3))") "Frequency", hb_count(1), hb_count(2), hb_count(3), hb_count(4)
	write(9,"(A12,4(f12.3))") "Percent", hb_count(1)*100/sum(hb_count), hb_count(2)*100/sum(hb_count),&
	hb_count(3)*100/sum(hb_count), hb_count(4)*100/sum(hb_count)
	
	close(8); close(9)

end program hbd

! Work subroutine for hydrogen bond counting in a given timestep
subroutine ts_hbd(natoms, labels, types, hbdist, box)
	implicit none
	
	! Variables comming from outside
	integer :: natoms, types
	real :: box(3), hbdist(types+1,5)
	character :: labels(types,2)*8
	
	! Subroutine Variables
	real :: rhb=2.67, coord(3), jHO(3), jOH(3), distHO, distOH
	logical :: box_corr(2,3)
	character :: label*8
	integer :: i, j, k, n, natoms, nmol(types), hbatoms, code, dim=3
	
	! Allocatable
	real :: hbmol(:,:), temp(:,:), iHO(:,:), iOH(:,:)
	
	nmol = 1
	hbatoms = 0
	allocate(hbmol(0,7))
	
	! Collect atoms of each molecule
	do k=1,natoms
		read(8,"(A8)") label
		do i=1,types
			
			if(label==labels(i,1)) then
				
				! Increment hbmol and nmol
				allocate(temp(size(hbmol,1)+1, size(hbmol,2)))
				temp(1:size(hbmol,1),:) = hbmol
				call move_alloc(temp, hbmol)
				
				nmol(i) = nmol(i) + 1
				
				! Read coords
				read(8,*) coord(1), coord(2), coord(3)
				
				! Put molecule type and coords in hbmol
				hbmol(size(hbmol,1),1:4) = (/ i, coord /)
				
				do
					read(8,"(A8)") label
					if(label==labels(i,2)) then
				
						! Read coords
						read(8,*) coord(1), coord(2), coord(3)
				
						! Put molecule type and coords in hbmol
						hbmol(size(hbmol,1),5:7) = coord
				
						exit ! go to next i
					else
						read(8,*)
					end if
				end do
			else
				read(8,*)
			end if
		end do
	end do
	
	! Count hydrongen bondings
	! Codes: 0 - none; 1 - hydrogen; 2 - oxygen; 3 - both;
	do i=1,sum(nmol)-1
		code = 0
		
		! Boundary conditions check
		! Determine how many projections to compare
		do k=1,size(box)
			if(abs(hbmol(i,1+k)) > box(k)/2-rhb ) then
				box_corr(1,k) = .true.
			else
				box_corr(1,k) = .false.
			end if

			if(abs(hbmol(i,4+k)) > box(k)/2-rhb ) then
				box_corr(2,k) = .true.
			else
				box_corr(2,k) = .false.
			end if
		end do
		
		allocate(iHO(2**count(box_corr(1,:)), dim), iOH(2**count(box_corr(2,:)), dim))
		
		do k=1,size(iHO,1) ! For iHO
			if(k==1) then ! Original
				iHO(k,1:3) = hbmol(i,2:4)
			else if(box_corr(1,1)) then ! X shift
				iHO(k,2:3) = hbmol(i,3:4)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbmol(i,2)-box(1)
				else
					iHO(k,1) = hbmol(i,2)+box(1)
				end if
			else if(box_corr(1,2)) then ! Y shift
				iHO(k,1) = hbmol(i,2)
				iHO(k,3) = hbmol(i,4)
				if(iHO(1,2) > 0) then
					iHO(k,2) = hbmol(i,3)-box(2)
				else
					iHO(k,2) = hbmol(i,3)+box(2)
				end if
				
			else if(box_corr(1,3)) then ! Z shift
				iHO(k,1:2) = hbmol(i,2:3)
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbmol(i,4)-box(3)
				else
					iHO(k,3) = hbmol(i,4)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2)) then ! XY shift
				iHO(k,3) = hbmol(i,4)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbmol(i,2)-box(1)
				else
					iHO(k,1) = hbmol(i,2)+box(1)
				end if
				if(iHO(1,2) > 0) then
					iHO(k,2) = hbmol(i,3)-box(2)
				else
					iHO(k,2) = hbmol(i,3)+box(2)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,3)) then ! XZ shift
				iHO(k,2) = hbmol(i,3)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbmol(i,2)-box(1)
				else
					iHO(k,1) = hbmol(i,2)+box(1)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbmol(i,4)-box(3)
				else
					iHO(k,3) = hbmol(i,4)+box(3)
				end if
				
			else if(box_corr(1,2) .and. box_corr(1,3)) then ! YZ shift
				iHO(k,1) = hbmol(i,2)
				if(iHO(1,2) > 0) then
					iHO(k,2) = hbmol(i,3)-box(2)
				else
					iHO(k,2) = hbmol(i,3)+box(2)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbmol(i,4)-box(3)
				else
					iHO(k,3) = hbmol(i,4)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2) .and. box_corr(1,3)) then ! XYZ shift
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbmol(i,2)-box(1)
				else
					iHO(k,1) = hbmol(i,2)+box(1)
				end if
				if(iHO(1,2) > 0) then
					iHO(k,2) = hbmol(i,3)-box(2)
				else
					iHO(k,2) = hbmol(i,3)+box(2)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbmol(i,4)-box(3)
				else
					iHO(k,3) = hbmol(i,4)+box(3)
				end if
			end if
		end do
		
		do k=1,size(iOH,1) ! For iOH
			if(k==1) then ! Original
				iOH(k,1:3) = hbmol(i,5:7)
			else if(box_corr(1,1)) then ! X shift
				iOH(k,2:3) = hbmol(i,6:7)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbmol(i,5)-box(1)
				else
					iOH(k,1) = hbmol(i,5)+box(1)
				end if
			else if(box_corr(1,2)) then ! Y shift
				iOH(k,1) = hbmol(i,5)
				iOH(k,3) = hbmol(i,7)
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbmol(i,6)-box(2)
				else
					iOH(k,2) = hbmol(i,6)+box(2)
				end if
				
			else if(box_corr(1,3)) then ! Z shift
				iOH(k,1:2) = hbmol(i,5:6)
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbmol(i,7)-box(3)
				else
					iOH(k,3) = hbmol(i,7)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2)) then ! XY shift
				iOH(k,3) = hbmol(i,7)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbmol(i,5)-box(1)
				else
					iOH(k,1) = hbmol(i,5)+box(1)
				end if
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbmol(i,6)-box(2)
				else
					iOH(k,2) = hbmol(i,6)+box(2)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,3)) then ! XZ shift
				iOH(k,2) = hbmol(i,6)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbmol(i,5)-box(1)
				else
					iOH(k,1) = hbmol(i,5)+box(1)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbmol(i,7)-box(3)
				else
					iOH(k,3) = hbmol(i,7)+box(3)
				end if
				
			else if(box_corr(1,2) .and. box_corr(1,3)) then ! YZ shift
				iOH(k,1) = hbmol(i,5)
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbmol(i,6)-box(2)
				else
					iOH(k,2) = hbmol(i,6)+box(2)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbmol(i,7)-box(3)
				else
					iOH(k,3) = hbmol(i,7)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2) .and. box_corr(1,3)) then ! XYZ shift
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbmol(i,5)-box(1)
				else
					iOH(k,1) = hbmol(i,5)+box(1)
				end if
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbmol(i,6)-box(2)
				else
					iOH(k,2) = hbmol(i,6)+box(2)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbmol(i,7)-box(3)
				else
					iOH(k,3) = hbmol(i,7)+box(3)
				end if
			end if
		end do
		write(*,*) size(iHO,1), size(iOH,1)
		
		! Determine, if any, what are the hydrogen bonds
		do j=i+1,nmol
			jHO = hbmol(j,4:6)
			jOH = hbmol(j,7:9)
			
			! Check H···O
			if(hbmol(j,2)==0) then
				do k=1,size(iOH,1)
					distHO = sqrt( sum( (jHO-iOH(k,1:3))**2 ) )
					if(distHO < rhb) then
						hbmol(j,2) = 1
						hbmol(i,3) = 1
					end if
				end do
			end if
			
			! Check O···H
			if(hbmol(j,3)==0) then
				do k=1,size(iHO,1)
					distOH = sqrt( sum( (jOH-iHO(k,1:3))**2 ) )
					if(distOH < rhb) then
						hbmol(j,3) = 1
						hbmol(i,2) = 1
					end if
				end do
			end if
		end do
		
		deallocate(iHO, iOH)
		
		! Add the molecule to the statistics
		code = hbmol(i,2) + hbmol(i,3)*2 + 1
		hb_count(code) = hb_count(code) + 1
	end do
end subroutine ts_hbd
			
			
		