! hbd
!
! Created on: Jun 2, 2014
! Author: Pedro Duarte and Pedro Morgado
! License: MIT License
! Description: Determine the distribution of hydrogen bonds in mixtures – both O and H participate
! in hydrogen bonds; only H; only O; and none (the molecule doesn't make any hydrogen bonds).

program hbd

	implicit none
	! Variable declaration
	integer, parameter :: molecules=500
	integer :: status=0, index, natoms
	character :: label*8, skip, labels(2)*8
	real :: hb_count(4), ts=0.0, box(3)
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-hbd", status="replace", action="write") ! File to write output
	
	! Removes the first two header lines
	read(8,*)
	read(8,*)
	
	! Get labels from command line
	labels = (/ "HOT", "OHT" /)
	hb_count = 0
	
	natoms = 4500
	do	
		read(8,*,iostat=status) label, index
		
		if(status==-1) then
			EXIT
		end if
		
		if(label=="timestep") then
			ts = ts + 1
			write(*,*) label, natoms, hb_count, labels
			
			read(8,*) box(1)
			read(8,*) skip, box(2)
			read(8,*) skip, skip, box(3)
			
			call ts_hbd(natoms, labels, hb_count, molecules, box)
		end if
	end do
	
	! Average frequencies
	hb_count = hb_count/ts
	
	write(9,"(5(A12))") "HBonds", "None", "Hydrogen", "Oxygen", "Both" 
	write(9,"(A12,4(f12.3))") "Frequency ", hb_count(1), hb_count(2), hb_count(3), hb_count(4)
	write(9,"(A12,4(f12.3))") "Percent ", hb_count(1)*100/sum(hb_count), hb_count(2)*100/sum(hb_count),&
	hb_count(3)*100/sum(hb_count), hb_count(4)*100/sum(hb_count)
	
	close(8); close(9)

end program hbd

! Work subroutine for hydrogen bond counting in a given timestep
subroutine ts_hbd(natoms, labels, hb_count, molecules, box)
	implicit none
	integer :: molecules
	real :: rhb=2.67, hbarray(molecules,9), coord(3), jHO(3), jOH(3), distHO, &
	distOH, hb_count(4), box(3)
	logical :: box_corr(2,3)
	real, allocatable :: iHO(:,:), iOH(:,:)
	character :: labels(2)*8, alabel*8
	integer :: i, j, k, natoms, nmol, hbatoms, status, code, dim=3
	
	nmol = 1
	hbatoms = 0
	
	! Collect atoms of each molecule
	do k=1,natoms
		read(8,"(A8)",iostat=status) alabel
		
		if(alabel == "HOT") then
			read(8,*) coord(1), coord(2), coord(3)
			hbatoms = hbatoms + 1
			
			hbarray(nmol,1:3) = (/ 1, 0, 0 /)
			hbarray(nmol,4:6) = coord
		
		else if(alabel == "OHT") then
			read(8,*) coord(1), coord(2), coord(3)
			hbatoms = hbatoms + 1
			
			hbarray(nmol,7:9) = coord
		
		else
			read(8,*)
		end if
		
		if(hbatoms == size(labels)) then
			hbatoms = 0
			nmol = nmol + 1
		end if
	end do
	
	! Determine if they have hydrogen bonding and where
	! Codes: 0 - none; 1 - hydrogen; 2 - oxygen; 3 - both;
	do i=1,nmol-1
		code = 0
		
		! Boundary conditions check
		! Determine how many projections to compare
		do k=1,size(box)
			if(abs(hbarray(i,3+k)) > box(k)/2-rhb ) then
				box_corr(1,k) = .true.
			else
				box_corr(1,k) = .false.
			end if

			if(abs(hbarray(i,6+k)) > box(k)/2-rhb ) then
				box_corr(2,k) = .true.
			else
				box_corr(2,k) = .false.
			end if
		end do
		
		allocate(iHO(2**count(box_corr(1,:)), dim), iOH(2**count(box_corr(2,:)), dim))
		
		do k=1,size(iHO,1) ! For iHO
			if(k==1) then ! Original
				iOH(k,1:3) = hbarray(i,4:6)
			else if(box_corr(1,1)) then ! X shift
				iHO(k,2:3) = hbarray(i,5:6)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbarray(i,4)-box(1)
				else
					iHO(k,1) = hbarray(i,4)+box(1)
				end if
			else if(box_corr(1,2)) then ! Y shift
				iHO(k,1) = hbarray(i,4)
				iHO(k,3) = hbarray(i,6)
				if(iHO(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				
			else if(box_corr(1,3)) then ! Z shift
				iHO(k,1:2) = hbarray(i,4:5)
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbarray(i,6)-box(3)
				else
					iHO(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2)) then ! XY shift
				iHO(k,3) = hbarray(i,6)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbarray(i,4)-box(1)
				else
					iHO(k,1) = hbarray(i,4)+box(1)
				end if
				if(iHO(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,3)) then ! XZ shift
				iHO(k,2) = hbarray(i,5)
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbarray(i,4)-box(1)
				else
					iHO(k,1) = hbarray(i,4)+box(1)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbarray(i,6)-box(3)
				else
					iHO(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,2) .and. box_corr(1,3)) then ! YZ shift
				iHO(k,1) = hbarray(i,4)
				if(iHO(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbarray(i,6)-box(3)
				else
					iHO(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2) .and. box_corr(1,3)) then ! XYZ shift
				if(iHO(1,1) > 0) then
					iHO(k,1) = hbarray(i,4)-box(1)
				else
					iHO(k,1) = hbarray(i,4)+box(1)
				end if
				if(iHO(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				if(iHO(1,3) > 0) then
					iHO(k,3) = hbarray(i,6)-box(3)
				else
					iHO(k,3) = hbarray(i,6)+box(3)
				end if
			end if
		end do
		
		do k=1,size(iOH,1) ! For iOH
			if(k==1) then ! Original
				iOH(k,1:3) = hbarray(i,4:6)
			else if(box_corr(1,1)) then ! X shift
				iOH(k,2:3) = hbarray(i,5:6)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbarray(i,4)-box(1)
				else
					iOH(k,1) = hbarray(i,4)+box(1)
				end if
			else if(box_corr(1,2)) then ! Y shift
				iOH(k,1) = hbarray(i,4)
				iOH(k,3) = hbarray(i,6)
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				
			else if(box_corr(1,3)) then ! Z shift
				iOH(k,1:2) = hbarray(i,4:5)
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbarray(i,6)-box(3)
				else
					iOH(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2)) then ! XY shift
				iOH(k,3) = hbarray(i,6)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbarray(i,4)-box(1)
				else
					iOH(k,1) = hbarray(i,4)+box(1)
				end if
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,3)) then ! XZ shift
				iOH(k,2) = hbarray(i,5)
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbarray(i,4)-box(1)
				else
					iOH(k,1) = hbarray(i,4)+box(1)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbarray(i,6)-box(3)
				else
					iOH(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,2) .and. box_corr(1,3)) then ! YZ shift
				iOH(k,1) = hbarray(i,4)
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbarray(i,6)-box(3)
				else
					iOH(k,3) = hbarray(i,6)+box(3)
				end if
				
			else if(box_corr(1,1) .and. box_corr(1,2) .and. box_corr(1,3)) then ! XYZ shift
				if(iOH(1,1) > 0) then
					iOH(k,1) = hbarray(i,4)-box(1)
				else
					iOH(k,1) = hbarray(i,4)+box(1)
				end if
				if(iOH(1,2) > 0) then
					iOH(k,2) = hbarray(i,5)-box(2)
				else
					iOH(k,2) = hbarray(i,5)+box(2)
				end if
				if(iOH(1,3) > 0) then
					iOH(k,3) = hbarray(i,6)-box(3)
				else
					iOH(k,3) = hbarray(i,6)+box(3)
				end if
			end if
		end do
		write(*,*) size(iHO,1), size(iOH,1)
		! Determine, if any, what are the hydrogen bonds
		do j=i+1,nmol
			jHO = hbarray(j,4:6)
			jOH = hbarray(j,7:9)
			
			! Check H···O
			if(hbarray(j,2)==0) then
				do k=1,size(iOH,1)
					distHO = sqrt( sum( (jHO-iOH(k,1:3))**2 ) )
					if(distHO < rhb) then
						hbarray(j,2) = 1
						hbarray(i,3) = 1
					end if
				end do
			end if
			
			! Check O···H
			if(hbarray(j,3)==0) then
				do k=1,size(iHO,1)
					distOH = sqrt( sum( (jOH-iHO(k,1:3))**2 ) )
					if(distOH < rhb) then
						hbarray(j,3) = 1
						hbarray(i,2) = 1
					end if
				end do
			end if
		end do
		
		deallocate(iHO, iOH)
		
		! Add the molecule to the statistics
		code = hbarray(i,2) + hbarray(i,3)*2 + 1
		hb_count(code) = hb_count(code) + 1
	end do
end subroutine ts_hbd
			
			
		