! surfangle
!
! Created on: Apr 28, 2014
! Author: Pedro Duarte and Miguel Teixeira
! License: MIT License
! Version: 0.3.1
! Description: Calculate the angle from the plane on the base of the 
! molecule, to a mean segment representing the molecule. Uses DL_POLY
! HISTORY files as input.

program surfangle

	implicit none
	! Variable declaration
	! Allocatable
	character, allocatable :: labels(:,:)*8
	real, allocatable :: global_angle(:), global_nmol(:), hist(:,:)
	! Allocated
	character :: arg*20, label*8, skip
	integer :: narg, index, status, frame, i, j, k, natoms, types
	real :: box(3), ts=0.
	
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
					write(*,*) "v0.3.1"
					stop
				case("--help","-h")
write(*,'(A10,(/))') "Surfangle"
write(*,*) "This utility calculates the angle to the surface of molecules forming a monolayer"
write(*,*) "on a surface. It also gives the orientation angle of the molecules. the program is designed"
write(*,*) "to work with any number of types of molecules on the monolayer. In adition, a histogram is"
write(*,*) "given of the frequency (in %) of the distribution of surface angles."
write(*,*) "Default run executes the program looking for only one molecule type, defined by the"
write(*,*) "terminal atoms OH (bottom) and CT (top)." 
write(*,'(A59,(/))') "To change the default labels use the option --labels (-l)."
write(*,'(A11,A69,A25,(/))') "--labels: ", &
"accepts an even number of labels in front; of each pair, the first is", " considered to be on top;"
write(*,'(A33,(/),A46,A82)') "Example: --labels OH CT OHA CTA;","Executes the program for two molecule types, ",&
"one with the base at OH and top at CT, the second with base at OHA and top at CTA."
stop
				case default
					write(*,*)"Option unknown: ", adjustl(arg)
					stop
			end select
		end do
	else
		allocate(labels(1,2))
		labels(1,1:2) = (/ "OH", "CT" /)
	end if
	
	types = size(labels,1)
	allocate(hist(types, 9), global_angle(types), global_nmol(types))
	hist = 0; global_angle = 0; global_nmol = 0;
	
	! File opening properties
	open(8, file="HISTORY", status="old", action="read") ! Open HISTORY
	open(9, file="HISTORY-angle", status="replace", action="write") ! Create and open output for average angles
	open(10, file="HISTORY-angles", status="replace", action="write") ! Create and open output for all angles
	open(11, file="HISTORY-hist", status="replace", action="write") ! Create and open output for histograms
	open(12, file="HISTORY-cris", status="replace", action="write") ! Create and open output for cristalinity
	
	! Files HEADER's
	write(9,*) "Timestep average angles to surface (Degrees). Average over all simulation at the end of file."
	write(10,*) "Angles to surface (surf) and orientation angles (ori). Angles are all in Degrees."
	write(11,*) "Angles distribution frequency (%)."
	write(11,'(A80,(/),A12,A12,9(I12))') "Angle interval","timestep","molecule",10,20,30,40,50,60,70,80,90
	write(12,'(A12,2X,A8)') "timestep", "cris"
	
	
	! Discard HEADER
	read(8,*); read(8,*) skip, skip, natoms;
	
	do
		read(8,*,iostat=status) label, index
		
		if(status .eq. -1) then
			EXIT
		end if
	
		if(label=="timestep") then
			ts = ts + 1
			read(8,*) box(1)
			read(8,*) skip, box(2)
			read(8,*) skip, skip, box(3)
			
			write(10,'(9(A12),(/),I12)') "timestep", "surf", "ori", "u(x)", "u(y)", "u(z)",&
			"v(x)", "v(y)", "v(z)", index
			
			call calc_angles(index, natoms, labels, types, box, global_nmol, global_angle, hist)
		end if
	end do

	! Write Average Angle
	write(9,'(A12,2(f12.4))') "Average", global_angle/global_nmol
	
	! Write Average Histogram
	do k=1,types
		hist(k,:) = hist(k,:)/global_nmol(k)*100
	end do
	write(11,'(A12,(/),12X,A12,9(f12.3))') "Average", trim(labels(1,1))//"-"//trim(labels(1,2)), hist(1,:)
	do i=2,types
		write(11,'(12X,A12,9(f12.3))') trim(labels(i,1))//"-"//trim(labels(i,2)), hist(i,:)
	end do
	
	close(8); close(9); close(10); close(11); close(12);

end program surfangle

subroutine calc_angles(index, natoms, labels, types, box, global_nmol, global_angle, hist)
	implicit none
	
	! Variables comming from outside
	integer :: types, index, natoms
	character :: labels(types,2)*8
	real :: box(3), global_angle(types), global_nmol(types), hist(types,9)
	
	! Subroutine variables
	! Arrays
	real, allocatable :: mol(:), temp(:), rdf(:)
	integer :: frame_hist(types,9)
	real :: u(3), v(3), ex(3)=[1,0,0], ez(3)=[0,0,1], nmol(types), angles(types), coord(types*2,3),&
	cm(2)=[0,0]
	logical :: check(types,2)
	! Single Value
	integer :: i, j, k, n, comb
	character :: label*8, skip
	real :: angle, ori_angle=0.0, tam, cris, cm_d
	! Parameters
	real, parameter :: pi = 4 * atan(1.0)
	
	! Initialize variables
	check = .false.
	angles = 0
	nmol = 0
	frame_hist = frame_hist*0
	cris = 0;
	allocate(mol(0)) ! start with empty array
	
	do n=1,natoms
		read(8,*) label
			
		outer: do i=1,types+1
			if(i==types+1) then
				read(8,*)
				
				exit
			end if
			do j=1,2
				if(label==labels(i,j)) then
					if(i .le. j) then
						read(8,'((2X,ES11.4E3),2(1X,ES11.4E3))') coord(i*j,1), coord(i*j,2), coord(i*j,3)
					else
						read(8,'((2X,ES11.4E3),2(1X,ES11.4E3))') coord(i+int(ceiling(real(i)/2)),1),&
						coord(i+int(ceiling(real(i)/2)),2), coord(i+int(ceiling(real(i)/2)),3)
					end if
					check(i,j) = .true.
					exit outer
				end if
			end do
		end do outer
		
		do k=1,types
			if(count(check(k,1:2))==2) then
				check(k,1:2) = .false.
				
				! Increment mol array
				allocate(temp(size(mol)+4))
				temp(size(mol)+1:size(temp)) = 0
				temp(1:size(mol)) = mol
				call move_alloc(temp, mol)
				! Increment rdf array
				allocate(temp(size(rdf)+3))
				temp(size(rdf)+1:size(temp)) = 0
				temp(1:size(rdf)) = rdf
				call move_alloc(temp, rdf)
				
				! Boundary conditions correction
				if(abs(coord(k,1)-coord(k+1,1)) > (box(1)/2)) then
					if( coord(k+1,1) < 0 ) then
						coord(k+1,1) = coord(k+1,1) + box(1)
					else          
						coord(k+1,1) = coord(k+1,1) - box(1)
					end if
				end if
			
				if(abs(coord(k,2)-coord(k+1,2)) > (box(2)/2)) then
					if( coord(k+1,2) < 0 ) then
						coord(k+1,2) = coord(k+1,2) + box(2)
					else
						coord(k+1,2) = coord(k+1,2) - box(2)
					end if
				end if
				
				! Vectorize and normalize
				u = 0
				tam = 0
				do j=1,3
					u(j) = coord(k*2,j)-coord(k*2-1,j)
				end do
				
				! RDF
				rdf(size(rdf)-2) = k
				rdf(size(rdf)-1:size(rdf)) = u(1:2)/3
				
				tam = sqrt(sum(u*u))
				do j=1,3
					u(j) = u(j)/tam
				end do
				v = (/ u(1:2),0. /)
				
				! Cristlinity calculation
				mol(size(mol)-3) = k
				mol(size(mol)-2:size(mol)) = u
				do i=1,size(mol)-4,4
					cris = cris + dot_product(mol(i+1:i+3),u)
				end do
						 
				
				! Calculate angles
				ori_angle = acos(dot_product(v,ex))*360/2/pi
				
				if(u(2) .lt. 0) then
					ori_angle = ori_angle + 180
				end if
			
				angle = acos(dot_product(u,ez))*360/2/pi
				angles(k) = angles(k) + angle
				nmol(k) = nmol(k) + 1
			
				! Write to all angles
				write(10,'(A12,8(f12.4))') trim(labels(k,1))//"-"//trim(labels(k,2)),angle,ori_angle,u,v
			
				! Create histogram
				do i=1,9
					if(angle .lt. i*10) then
						frame_hist(k,i) = frame_hist(k,i) + 1
					end if
				end do
			end if
		end do
	end do
	
	! Determine rdf center
	do i=1,2
		do j=1,size(rdf),4
			cm(i) = cm(i) + rdf(i+j)
		end do
	end do
	cm = cm/(size(rdf)/3)
	
	do i=1,10
		do j=1,size(rdf),3
			cm_d = sqrt((rdf(j+1:j+2)-cm)**2)
			
			if(cm_d .gt. (i-1)*box(1)/10 .and. cm_d .lt. i*box(1)/10) then
				rdfs(i) = rdfs(i) + dot_product(mol)
	
	! Write frame averages
	write(9,'(I12,2(f12.4))') index, angles/nmol

	! Histogramfor the frame
	hist = hist + frame_hist
	write(11,'(I12,A12,9(f12.3))') index, trim(labels(1,1))//"-"//trim(labels(1,2)), frame_hist(1,:)*100/nmol(1)
	do i=2,types
		write(11,'(12X,A12,9(f12.3))') trim(labels(i,1))//"-"//trim(labels(i,2)), frame_hist(i,:)*100/nmol(i)
	end do

	! Global Average Angle
	global_angle = global_angle + angles
	global_nmol = global_nmol + nmol
	
	! Write Cristalinity
	write(12,"(I12,2X)",advance="no") index;write(12,*) cris/comb(size(mol)/4);
	deallocate(mol)
end subroutine calc_angles

! Calculate the number of combinations
function comb(n) result(r)
	integer, intent(in) :: n
	integer :: r, i
	
	r = 0
	do i=1,n
		r = r+(i-1)
	end do
	
end function comb