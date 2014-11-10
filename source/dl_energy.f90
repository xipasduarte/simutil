! dl_energy
!
! Created on: Jul 30, 2014
! Author: Pedro Duarte
! License: MIT License
! Version: 0.0.1
! Description: Extract data from DL_POLY's OUTPUT or STATIS file.

program dl_energy
	implicit none
	! Variable declaration
	! Allocatable
	real, allocatable :: values(:), output(:)
	integer, allocatable :: properties(:), temp_i(:)
	! Allocated
	character :: arg*20, var*2, label_list(30)*6, keys*70, srcfile*80="STATIS", outfile*80="STATIS-out", line*80
	integer :: status, narg, i, j, k, ts_tot=0, ts, arr_size, prop_size
	real :: time
	
	! Handle command line arguments
	! Check if any arguments are found
	narg = command_argument_count()
	! Loop over the arguments
	if(narg>0)then
		! loop across options
		do i=1,narg,2
			call get_command_argument(i,arg)
			select case(adjustl(arg))
				case("--srcfile", "-s")
					call get_command_argument(i+1,arg)
					srcfile=arg
				case("--outfile", "-o")
					call get_command_argument(i+1,arg)
					outfile=arg
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
	
	label_list=(/ "engcns", "temp  ", "engcfg", "engsrp", "engcpe", "engbnd", "engang", "engdih", "engtet", &
	"enthal", "tmprot", "vir   ", "virsrp", "vircpe", "virbnd", "birang", "vircon", "virtet", "volume", "tmpshl", &
	"engshl", "virshl", "alpha ", "beta  ", "gamma ", "virpmf", "press ", "amsd  ", "stress", "cell  " /)
	
	do i=1,size(label_list, 1)
		if(mod(i,5)==0) then
			write(*,"(i2,a3,a6,2X)", advance="yes") i, " - ", label_list(i)
		else
			write(*,"(i2,a3,a6,2X)", advance="no") i, " - ", label_list(i)
		end if
	end do
	write(*,"(/,a80)") adjustl("Select values to extract (comma separated).")
	read(*,"(a70)") keys
	
	! Abort program if it has no keys
	if(len_trim(keys)==0) then
		write(*,*) "No properties were selected."
		call abort
	end if
	
	! Convert keys string into an array
	allocate(properties(1)); properties(:) = 0;
	
	do i=1,len_trim(keys)
		prop_size=size(properties,1)
		if(keys(i:i)/=",") then
			if(properties(size(properties,1))==0) then
				read(keys(i:i), "(i2)") properties(prop_size)
			else
				var = char(properties(prop_size))//keys(i:i)
				read(var, "(i2)") properties(prop_size)
			end if
		else ! Increment properties array
			allocate(temp_i(prop_size+1))
			temp_i(1:size(properties,1)) = properties(:)
			temp_i(size(temp_i,1)) = 0
			call move_alloc(temp_i, properties)
		end if
	end do
	
	! Allocate output string with required size
	allocate(output(size(properties,1)))
	
	! Open files
	open(1, file=srcfile, status="old", action="read") ! Open STATIS
	open(2, file=outfile, status="replace", action="write") ! Open output file
	
	! Transcribe header
	read(1,"(a80)") line; write(2,*) adjustl(line);
	read(1,"(a80)") line; write(2,*) adjustl(line);
	! Write labels
	write(2,"(a10,a14)", advance="no") adjustl("timestep"), adjustl("time")
	do i=1,size(properties,1)
		if(i==size(properties,1)) then
			write(2,"(a14)", advance="yes") adjustl(label_list(properties(i)))
		else
			write(2,"(a14)", advance="no") adjustl(label_list(properties(i)))
		end if
	end do
	
	do
		! Read timestep header
		read(1,"(i10,e14.6,i10)", iostat=status) ts, time, arr_size
		
		! Exit if there are no more timesteps
		if(status==-1) then
			EXIT
		end if
		
		! Write timestep information
		write(2,"(i10,e14.6)", advance="no") ts, time
		
		! Get all values into values() array
		allocate(values(arr_size))
		do i=1,arr_size
			if(mod(i,5)==0 .or. i==arr_size) then
				read(1,"(e14.6)", advance="yes") values(i)
			else
				read(1,"(e14.6)", advance="no") values(i)
			end if
		end do
		
		! Select requested values and record in output file and array
		do i=1,size(properties,1)
			if(i==size(properties,1)) then
				write(2,"(e14.6)", advance="yes") values(properties(i))
			else
				write(2,"(e14.6)", advance="no") values(properties(i))
			end if
			
			output(i) = output(i) + values(properties(i))
		end do
		deallocate(values)
	end do
					
end program dl_energy