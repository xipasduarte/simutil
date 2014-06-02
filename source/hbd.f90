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
	integer :: n, v
	parameter ( n=600, v=8 )
	real, dimension(n*v) :: hbarray
	integer :: nmol, status=0, i, j, k, index
	character :: label*8, skip*80, arg_name*20
	real :: xO=0, yO=0, zO=0, xH=0, yH=0, zH=0, rhb=2.67, distOH, distHO, oh_hb=0.0, o_hb=0.0, h_hb=0.0,&
	no_hb=0.0, ts=0.0
	
	open(8, file="HISTORY", status="old", action="read")
	open(9, file="HISTORY-hbd", status="replace", action="write") ! File to write output
	
	! Removes the first two header lines
	read(8,*)
	read(8,*)
	
	! Get labels from command line
	
	
	do	
		read(8,*,iostat=status) label, index
		
		if(status==-1) then
			EXIT
		end if
		
		if(label=="timestep") then ! Discard timestep and box information
			read(8,*)
			read(8,*)
			read(8,*)
			
			if(ts > 0) then
				do i=1,nmol
					if(hbarray((i-1)*v+1) == 1 .and. hbarray((i-1)*v+2) == 1) then
						oh_hb = oh_hb + 1
					else if(hbarray((i-1)*v+1) == 1) then
						h_hb = h_hb + 1
					else if(hbarray((i-1)*v+2) == 1) then
						o_hb = o_hb + 1
					else
						no_hb = no_hb + 1
					end if
				end do
			end if
			
			ts = ts+1
			nmol = 0
		
		else if(label=="HOT") then ! Get hydrogen coord
			read(8,*) xH, yH, zH
			
		else if(label=="OHT") then ! Get oxygen coord
			read(8,*) xO, yO, zO
			
		else
			read(8,*,iostat=status)
		end if
		
		if(xH.ne.0 .and. yH.ne.0 .and. zH.ne.0 .and. xO.ne.0 .and. yO.ne.0 .and. zO.ne.0) then
			nmol = nmol + 1

			hbarray(nmol*v+1) = 0
			hbarray(nmol*v+2) = 0
			hbarray(nmol*v+3) = xH
			hbarray(nmol*v+4) = yH
			hbarray(nmol*v+5) = zH
			hbarray(nmol*v+6) = xO
			hbarray(nmol*v+7) = yO
			hbarray(nmol*v+8) = zO
			
			do i=1, nmol
				if(i == nmol) then
					
					! Reset coord
					xH = 0; yH = 0; zH = 0; xO = 0; yO = 0; zO = 0;
					
					EXIT
				end if
				
				if(hbarray((i-1)*v+1)==0) then
					distHO = sqrt( (hbarray((i-1)*v+3)-xO)**2 + (hbarray((i-1)*v+4)-yO)**2 + (hbarray((i-1)*v+5)-zO)**2 )
					if(distHO < rhb) then
						hbarray((i-1)*v+1) = 1
						hbarray(nmol*v+2) = 1
					end if
				end if
				
				if(hbarray((i-1)*v+2)==0) then
					distOH = sqrt( (hbarray((i-1)*v+6)-xH)**2 + (hbarray((i-1)*v+7)-yH)**2 + (hbarray((i-1)*v+8)-zH)**2 )
					if(distOH < rhb) then
						hbarray((i-1)*v+2) = 1
						hbarray(nmol*v+1) = 1
					end if
				end if
			end do
		end if
		
		if(status==-1) then
			EXIT
		end if
	end do
	
	! Final count
	do i=1,nmol
		if(hbarray((i-1)*v+1) == 1 .and. hbarray((i-1)*v+2) == 1) then
			oh_hb = oh_hb + 1
		else if(hbarray((i-1)*v+1) == 1) then
			h_hb = h_hb + 1
		else if(hbarray((i-1)*v+2) == 1) then
			o_hb = o_hb + 1
		else
			no_hb = no_hb + 1
		end if
	end do
	
	write(9,"(5(A12))") "HBonds", "Both", "Hydrogen", "Oxygen", "None" 
	write(9,"(A12,4(f12.3))") "Frequency ", oh_hb/ts, h_hb/ts, o_hb/ts, no_hb/ts
	
	close(8); close(9)

end program hbd
