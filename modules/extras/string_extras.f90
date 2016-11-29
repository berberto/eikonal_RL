
module string_extras
    use constants
    implicit none

contains

    subroutine remove_blanks(string)
	    character(len=*) :: string
	    integer :: stringLen 
	    integer :: last, actual
	
	    stringLen = len (string)
	    last = 1
	    actual = 1
	
	    do while (actual < stringLen)
	        if (string(last:last) == ' ') then
	            actual = actual + 1
	            string(last:last) = string(actual:actual)
	            string(actual:actual) = ' '
	        else
	            last = last + 1
	            if (actual < last) &
	                actual = last
	        endif
	    end do
    end subroutine remove_blanks


	! Replaces white spaces with given character
    subroutine replace_blanks(string,ctr)
	    character(len=*) :: string
	    character(len=1) :: ctr		! character to replace to blanks
	    integer :: stringlen 
	    integer :: i, j, first, last

	    stringlen = len(string)	    

		! An index runs untill a white space is encountered;
		!	the position of the first white character is saved as "first".
		! Then index runs untill a characther non-white-space is found:
		!	the position of the last white space is saved as "last",
		!	and the portion of white spaces between "first" and "last" is replaced by a string of 'ctr'.
		! If "last" is the last character of the string, then do not replace and leave white.
		i = 1
	    do while (i < stringlen)
	        if (string(i:i) == ' ') then
				! find first and last white spaces
				first = i
				do while ((string(i:i) == ' ') .and. (i<stringlen))
					last = i
					i = i+1
				end do
				! replace white spaces with 'ctr'
				if (i .ne. stringlen) then
					do j=first, last
						string(j:j) = ctr
					end do
				else
					! if the white spaces are at the end of the string, do nothing
					exit
				end if
	        else
				i = i+1
			end if
	    end do
    end subroutine replace_blanks
    
end module string_extras
