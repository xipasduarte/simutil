! gro2cfg
!
! Created on: Mar 9, 2012
! Author: Pedro Duarte and Pedro Morgado
! License: MIT License
! Description: Convert a GROMACS .gro file into a DL Poly CONFIG file.

program gro2cfg

	implicit none
	! Variable declaration
	character(len=80) :: filename, outfilename, name, molecule, fmt, fmt_box
	character (len=132) :: skip
	character(len=8):: atom
	integer :: line, totlines, indexno
	real :: xbox, ybox, zbox, ax, ay, az, nullval
        nullval=0

	!Ask user for the path to the file to convert and output filename
	!print *, "Enter file path without the file (can be relative to this file location):"
	!read *, filepath
	filename = ""
    do while (filename == "")
        print *, "Enter filename (mustn't be null):"
        read *, filename
    end do
    print *, "Enter output filename:"
    read *, outfilename
	! Set filename to CONFIG if null input
	if (outfilename == "") filename="CONFIG"


	open(7, file=filename, status="old", action="read") ! To determine number of lines
	! Find number of lines of file and get box dimensions
	do line=1, 3
	    if (line == 1) then
	        read(7,*) name
	    else if (line == 2) then
	        read(7,*) totlines
	        totlines = totlines + 2 + 1 !Total lines is equal to all atoms plus 2 on top and 1 at the bottom
	    end if
	end do
	do line=3, totlines
	    if (line == totlines) then
            read(7,*) xbox, ybox, zbox
            xbox = xbox*10
            ybox = ybox*10
            zbox = zbox*10
        else
            read(7,*) skip
        end if
    end do

	close(7)
	open(8, file=filename, status="old", action="read") ! To do the actual conversion
	open(9, file=outfilename, status="replace", action="write") ! File to write output
	
	! Read positions and make CONFIG file
	fmt="(3(f20.12))"
    fmt_box="(3(f20.12))"
	!File HEADER
	write(9,'(A80)') name
	write(9,*) "        0","         1"
	write(9,fmt_box) xbox, nullval, nullval
        write(9,fmt_box) nullval, ybox, nullval
        write(9,fmt_box) nullval, nullval, zbox

	!File BODY
	do line=1, totlines
        if (line == 1 .or. line == 2 .or. line == totlines) then
            read(8,*) skip
        else
            read(8,*) molecule, atom, indexno, ax, ay, az
            ! Convert to Angstroms
            ax = ax*10
            ay = ay*10
            az = az*10

            if (ax<0) then
                ax = ax+xbox
            else if (ax>xbox) then
                ax = ax-xbox
            end if

	        if (ay<0) then
                ay = ay+ybox
            else if (ay>ybox) then
                ay = ay-ybox
            end if

            if (az<0) then 
                az = az+zbox
            else if (az>zbox) then
                az = az-zbox
            end if

            ax = ax-(xbox/2)
            ay = ay-(ybox/2)
            az = az-(zbox/2)
            write(9,'(A8,I10)') atom, indexno
            write(9,fmt) ax, ay, az
        end if
	end do

	close(8)
	close(9)

end program gro2cfg
