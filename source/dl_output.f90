! dl_energy
!
! Created on: Jul 30, 2014
! Author: Pedro Duarte
! License: MIT License
! Version: 0.0.1
! Description: Extract data from DL_POLY's OUTPUT or STATIS file.

program dl_output
	implicit none
	! Variable declaration
	! Allocatable
	real, allocatable :: values(:)
	integer, allocatable :: keys_arr(:), temp_i(:)
	! Allocated
	character :: arg*20, label_list(27)*12, keys*70, srcfile*80="OUTPUT", outfile*80="OUTPUT-dl_output", line*80,&
	label*9
	integer :: status, var, narg, i, j, k, ts_tot=0, ts, prop_size, stack=1, label_check, skiper=1
	logical :: col=.false.
	real :: mean(2,27)
	
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
				case("--col", "-c")
					call get_command_argument(i+1,arg)
					col=.true.
					keys=arg
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
	
	label_list=(/ "     eng_tot", "    temp_tot", "     eng_cfg", "     eng_vdw", "     eng_cou", "     eng_bnd",&
	"     eng_ang", "     eng_dih", "     eng_tet",	"      eng_pv", "    temp_rot", "    vir_conf",	"     vir_vdw",&
	"     vir_cou", "     vir_bnd", "     vir_ang",	"     vir_con", "     vir_tet", "      volume", "    temp_shl",&
	"     eng_shl", "     vir_shl", "       alpha", "        beta", "       gamma", "     vir_pmf", "       press" /)
	
	! Write properties to choose from on screen and get the user's selection
	if(.not. col) then
		do i=1,size(label_list, 1)
			if(mod(i,5)==0) then
				write(*,"(i2,a3,a12)", advance="yes") i, " - ", adjustl(label_list(i))
			else
				write(*,"(i2,a3,a12)", advance="no") i, " - ", adjustl(label_list(i))
			end if
		end do
		write(*,"(/,a42)") adjustl("Select values to extract (comma separated).")
		read(*,"(a70)") keys
	end if
	
	! Abort program if it has no keys
	if(len_trim(keys)==0) then
		write(*,*) "No properties were selected."
		call abort
	end if
	
	! Convert keys string into an array
	allocate(keys_arr(1)); keys_arr(:) = 0;
	do i=1,len_trim(keys)
		prop_size=size(keys_arr,1)
		if(keys(i:i)/=",") then
			if(keys_arr(prop_size)==0) then
				read(keys(i:i), "(i2)") keys_arr(prop_size)
			else
				read(keys(i:i),"(i2)") var
				keys_arr(prop_size) = keys_arr(prop_size)*10 + var
			end if
		else ! Increment keys_arr array
			allocate(temp_i(prop_size+1))
			temp_i(1:prop_size) = keys_arr(:)
			temp_i(size(temp_i,1)) = 0
			call move_alloc(temp_i, keys_arr)
		end if
	end do
	
	! Open files
	open(1, file=srcfile, status="old", action="read") ! Open OUTPUT
	open(2, file=outfile, status="replace", action="write") ! Open output file
	
	! Transcribe header & search for data stacking interval
	read(1,"(9(/),a80,13(/),30X,i12,8(/),30X,i12)") line, ts_tot, stack
	write(2,*) adjustl(line)
	
	! Write labels
	write(2,"(a9,2X)", advance="no") adjustl("step")
	do i=1,size(keys_arr,1)-1
			if(keys_arr(i) .lt. size(label_list,1)) then
				write(2,"(a12)", advance="no") adjustr(label_list(keys_arr(i)))
			else
				write(2,"(i12)", advance="no") keys_arr(i)
			end if
	end do
	if(keys_arr(size(keys_arr,1)) .lt. size(label_list,1)) then
		write(2,"(a12)", advance="yes") adjustr(label_list(keys_arr(size(keys_arr,1))))
	else
		write(2,"(i12)", advance="yes") keys_arr(size(keys_arr,1))
	end if
	
	! Search for print interval and tot_ts
	ts = stack
	
	! Process each timestep
	do while(ts_tot>=ts)
		! Search for " -------"
		read(1,"(a8)", iostat=status) label
		
		! Exit if there are no more timesteps
		if(status==-1) then
			EXIT
		end if
		
		if(label .eq. " -------") then
			if (skiper == 8) then
				skiper = 0; read(1,"(5/)")
			else
				! Read next line's first 8 characters and check to see if it is the time step
				read(1,"(i9)") label_check
			end if
			
			if (label_check == ts) then
				! Write and update ts and update skiper
				write(2,"(i9,2X)", advance="no") ts
				ts = ts + stack
				skiper = skiper + 1
				
				! Jump step values to averages
				read(1,"(2/)")
				
				allocate(values(size(label_list,1)))
				do i=1,3
					do j=1,9
						if(j==1) then
							read(1,"(9X,e12.4)", advance="no") values((i-1)*3+j)
						else
							read(1,"(e12.4)", advance="no") values((i-1)*3+j)
						end if
					end do
					read(1,*)
				end do
		
				! Select requested values and record in output file and array
				do i=1,size(keys_arr,1)
					if(i==size(keys_arr,1)) then
						write(2,"(e12.4)", advance="yes") values(keys_arr(i))
					else
						write(2,"(e12.4)", advance="no") values(keys_arr(i))
					end if
				end do
				deallocate(values)
				
			end if
		end if
	end do
	
	! Get Averages and R.M.S. fluctuations
	read(1,"(16/)") ! skip lines
	do k=1,2
		do i=1,3
			do j=1,9
				if(j==1) then
					read(1,"(9X,e12.4)", advance="no") mean(k,(i-1)*3+j)
				else
					read(1,"(e12.4)", advance="no") mean(k,(i-1)*3+j)
				end if
			end do
			read(1,*)
		end do
		read(1,*)
	end do

	! Write averages
	write(2,"(/,a9,2X)", advance="no") "Averages"
	do i=1,size(keys_arr,1)
		if(i==size(keys_arr,1)) then
			write(2,"(e12.4)", advance="yes") mean(1,keys_arr(i))
		else
			write(2,"(e12.4)", advance="no") mean(1,keys_arr(i))
		end if
	end do
	
	close(1); close(2)			
end program dl_output