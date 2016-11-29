module randomf

	use kinds, only: dp
	
	implicit none

contains

	subroutine seed_from_urandom(i)
	    implicit none
	    integer :: n,k
	    integer(4), intent(out), allocatable :: i(:)
	    real(dp) :: r

		call random_seed(size=n)
	    allocate(i(n))
	    open(89, file='/dev/urandom', access='stream', form='UNFORMATTED')
	    read(89) i
	    close(89)
	
	    call random_seed(put=i)

	end subroutine seed_from_urandom


	function random_integer (l, u)
		implicit none
		real :: r
		integer :: l, u
		integer :: random_integer

		call random_number(r)
		random_integer = floor(l + (u-l+1)*r)
		
	end function random_integer


end module randomf
