!
!	Module with constants and functions
!
module consts_funcs

	use kinds, only: dp
	use constants
	use mpi
	
	implicit none

!	include 'mpif.h'

	real(dp), parameter :: k = 1.				! curvature harmonic potential
	real(dp), parameter :: D = 1				! diffusion constant
	real(dp), parameter :: U = 2.				! speed
	real(dp), parameter :: Dx = .05				! lattice spacing
	real(dp), parameter :: Dt = .0001			! time discretization
	real(dp), parameter :: R = 10.				! square domain (-R,R)x(-R,R)
	integer, parameter :: L = ceiling(R/dx)		! lattice (-L,L)x(-L,L)
	integer, parameter :: T	= int(100./dt)		! maximum time
	integer, parameter :: Trlx = int(30./dt)	! time for equilibration 
	integer, parameter :: its = 5				! number of iterations

	real(dp), parameter :: delta = U*Dt/Dx
	real(dp), parameter :: sigma = 2.*D*Dt/(Dx**2.)

	private :: sigma, delta, R

	public :: Dx, Dt, k, L, T, Trlx, U, D, its

contains

	subroutine check_params ()
		if((sigma < delta).or.(sigma > .5)) then
			print *, "Check parameters, Dx and Dt"
			call exit(1)
		end if
	end subroutine check_params


	! policy function
	function policy (theta)
		real(dp) :: theta
		real(dp) :: policy(0:4)

		policy(0) = 1.- 2.*sigma
		policy(1) = .5*(sigma + delta*sin(theta))
		policy(2) = .5*(sigma + delta*cos(theta))
		policy(3) = .5*(sigma - delta*sin(theta))
		policy(4) = .5*(sigma - delta*cos(theta))
	end function policy


	! gradient of the log og the policy
	function grad_log_pi (theta)
		real(dp) :: theta
		real(dp) :: grad_log_pi(0:4)

		grad_log_pi(0) = 0.
		grad_log_pi(1) = delta*cos(theta)/(sigma + delta*sin(theta))
		grad_log_pi(2) = -delta*sin(theta)/(sigma + delta*cos(theta))
		grad_log_pi(3) = -delta*cos(theta)/(sigma - delta*sin(theta))
		grad_log_pi(4) = delta*sin(theta)/(sigma - delta*cos(theta))
		
	end function grad_log_pi


end module consts_funcs



!
!	MAIN PROGRAM
!
program eikonal2d

	use kinds, only: dp
	use constants, only: pi
	use string_extras
	
	use randomf
	use consts_funcs
	
	implicit none

	! MPI variables
	integer :: proc, Nproc
	integer :: ierr, error
	
	! command line inputs
	character(len=32) :: arg
	integer :: Na						! number of agents
	real(dp) :: g						! coupling strength
	real(dp) :: alpha, beta, eta		! learning rates
	real(dp) :: gvals(0:19)

	! random number
	integer, allocatable :: seed(:)
	real(dp) :: r

	! learning variables
	real(dp) :: p(0:4), elgb(0:4)			! policy and eligibility
	integer :: i, j, ii, jj, a, ts, it		! running indices
	integer :: x(2), x_new(2)
	integer, allocatable :: state(:,:), state_new(:,:) ! first index is the agent, second is the lattice position
	integer :: occupation(-L:L, -L:L), occupation_new(-L:L, -L:L)	! matrix with occupation numbers
	real(dp) :: rho_av(-L:L, -L:L)		! density averaged over several steps

	integer :: visits(-L:L, -L:L)		! counting of visits at a lattice point
	
	real(dp) :: v(-L:L, -L:L), v_new(-L:L, -L:L), v_av(-L:L,-L:L) ! table with the value of the state
	integer :: act(0:4, 2)	! first index is the action, second is lattice displacement
	real(dp) :: theta(-L:L, -L:L), theta_new(-L:L, -L:L) ! table of the parameter of the policy
	real(dp) :: drift_block(2)
	integer :: blocksize = 2
	real(dp) :: rew, rewbar, rewbar_new, avrew, delta
	real(dp) :: radial_rho(0:L), radial_val(0:L)

	real(dp) :: r1(2), r2(2)
	
	! input
	integer :: in_states=151, in_vt=152, states_error=153, vt_error=154
	character(len=62) :: in_form, in_name_states, in_name_vt

	! output
	integer, parameter :: Ntraj = 100
	integer :: out_std = 100
	integer :: out_traj(Ntraj) = 1000, out_cost=2001, out_theta=2002, &
			out_rho_value=2003, out_rho_av=2004, out_radial=2005, out_dev=2006
	character(len=64) :: out_form, out_dir, out_name, mk_dir

	do i=1,Ntraj
		out_traj(i) = out_traj(i) + i
	end do

	! check that policy is well defined
	call check_params

	!
	!	Parameters from command line
	! 
	if ((iargc() > 21) .or. (iargc() < 1)) then
		print *, "Input"
		print *, "1      -> number of agents"
		print *, "2(-21) -> array with values of the interaction strength"
		call exit (1)
	end if
	
	call getarg(1,arg)
	read(arg,*) Na

	do i=0, iargc()-2
		call getarg(i+2,arg)
		read(arg,*) gvals(i)
	end do

	!
	!	Initialize MPI environment
	!
call MPI_Init ( ierr )
	call MPI_Comm_size(MPI_COMM_WORLD, Nproc, error)
	call MPI_Comm_rank(MPI_COMM_WORLD, proc, error)

	if((Nproc > 20).and.(proc .eq. 0)) then
		print *, "Run with maximum 20 processors:"
		print *, "  mpirun -np <#proc> ./eikonal <#agents> <string of #proc values of g>"
		call exit(1)
	end if

	! # of process labels the value of g used in the simulation
	g = gvals(proc)
	alpha = .05
	beta = .8
	eta = .05

!	call  MPI_Finalize (ierr)
!	stop

	allocate(state(Na,2), state_new(Na,2))

	!
	!	Set the INPUT and OUTPUT files
	!
	call set_IO ()

	write(100+proc,fmt="(a, i2, a, f7.4, a)") "Process ", proc, ", g = ", g, " --> started..."
	
	! inizialize random number generator
	call seed_from_urandom(seed)

	! define the actions (displacements on the lattice)
	act(0,:) = (/ 0, 0 /)	! stay in place
	act(1,:) = (/ 0, 1 /)	! move UP
	act(2,:) = (/ 1, 0 /)	! move RIGHT
	act(3,:) = (/ 0, -1 /)	! move DOWN
	act(4,:) = (/ -1, 0 /)	! move LEFT

	! Initialization of occupation and visits matrices
	occupation_new = 0
	visits = 0

	! Initialization estimate of the stationary density
	rho_av = 0.
	v_av = 0.

	! Initialization VALUE function
	v = 0.	
	v_new = 0.
	rewbar = 0.
	rewbar_new = 0.
	
	!
	!	INITIALIZATION OF STATES AND POLICY
	!
	! --------------------------------------------------------------------------
	if ((states_error == 0).and.(vt_error == 0)) then
		! read walker positions from file
		do i=1, Na
			read(in_states,*) state(i,:)
		end do
		! read value and policy parameters from file
		read(in_vt,*) rewbar
		do i=-L,L
			do j=-L, L
				read(in_vt,*) v(i,j), theta_new(i,j)
			end do
		end do
	else
		! random initialization of the position of the walkers
		do i=1, Na
			state(i,:) = (/ random_integer(-L,L), random_integer(-L,L) /)
		end do
		! random initialization of the policy parameters
		do i=-L, L
			do j=-L, L
				call random_number(r)
				theta_new(i,j) = 2.*pi*r
			end do
		end do
	end if
	do i=1, Na
		occupation_new(state(i,1), state(i,2)) = occupation_new(state(i,1), state(i,2)) + 1
		visits(state(i,1), state(i,2)) = visits(state(i,1), state(i,2)) + 1
	end do
	close(in_states)
	close(in_vt)
	occupation = occupation_new


!	do i=-L, L
!		do j=-L, L
!			r1 = (/ -1.*i, -1.*j /)
!			r2 = (/ cos(theta_new(i,j)), sin(theta_new(i,j)) /)
!			r = acos(dot_product(r1,r2)/sqrt(sum(r1**2)))
!			write(100,*) i*dx, j*dx, r/pi
!		end do
!		write(100,*) ''
!	end do
!	stop

	!
	!	START THE LEARNING
	!
	! --------------------------------------------------------------------------
	! Take T steps for each iteration, for its iterations
	do it=1, its
		write(101+proc,fmt="(a, f7.4, a, i2)") "... g = ", g, ",  iteration ", it
		do ts=0, T
			! update theta field
			theta = theta_new
	
			!
			!	SNAPSHOTS every 10 time units, for the last iteration
			!
			! from time to time...
			if (mod(ts,int(10/dt))==0) then
				
				! .. print the drift vector field
				do i=-L, L-blocksize, blocksize
					do j=-L, L-blocksize, blocksize
						! write one direction vector per block (skipping the others)
!						drift_block = (/ cos(theta(i,j)), sin(theta(i,j)) /)
						! write the average over blocks of the direction vector
						drift_block = 0.
						do ii=0, blocksize-1
							do jj=0, blocksize-1
								drift_block(1) = drift_block(1) + cos(theta(i+ii, j+jj))
								drift_block(2) = drift_block(2) + sin(theta(i+ii, j+jj))
							end do
						end do
						drift_block = drift_block/sqrt(sum(drift_block**2.))

						r1 = (/ -1.*i, -1.*j /)
						r = acos(dot_product(r1,drift_block)/sqrt(sum(r1**2)))
						write(out_theta,*) i*dx, j*dx, (.5*blocksize)*dx*drift_block, r/pi, it
					end do
					write(out_theta,*) ''
				end do
				write(out_theta,*) ''
				write(out_theta,*) ''
	
				! .. the empirical density, the value ..
				do i=-L, L
					do j=-L, L
						write(out_rho_value,*) real(i*dx,4), real(j*dx,4), &
							real(1./(Na*dx**2)*occupation(i,j)), real(v(i,j),4), it
					end do
					write(out_rho_value,*) ''
				end do
				write(out_rho_value,*) ''
				write(out_rho_value,*) ''
	
				! .. and their radial averages
				radial_rho = radial_average(real(occupation, dp))/(Na*dx**2)
				radial_val = radial_average(v)
				do i=0, L
					write(out_radial,*) i*dx, radial_rho(i), radial_val(i), it
				end do
				write(out_radial,*) ''
				write(out_radial,*) ''
	
			end if

			!
			!	At each time step, move particles one at a time
			!	and update value and policy syncronously
			!
			avrew = 0
			do i=1, Na
				
				if ((mod(ts,int(.01/dt))==0).and.(i .le. Ntraj).and.(it==its)) then
					write(out_traj(i),*) ts*dt, dx*state(i,:)
				end if
				
				x = state(i,:)
				p = policy(theta(x(1), x(2)))
				elgb = grad_log_pi(theta(x(1), x(2)))
	
				! choose 'a' according to the policy
				call choose_action(a)
	
				! take action: move in the opposite direction if going outside the domain
				x_new(:) = x(:) + act(a,:)
				if ((abs(x_new(1)) .gt. L) .or. (abs(x_new(2)) .gt. L)) then
					x_new(:) = x(:) - act(a,:)
				end if
				
				rew = -(q(x_new) + collision(x_new))
	
				avrew = avrew + rew/Na
	
				delta = rew - rewbar + v(x_new(1),x_new(2)) - v(x(1), x(2))
	
				rewbar_new = rewbar_new + eta/(1.+(ts+1)**.6)*delta/Na
				
				v_new(x(1), x(2)) = v_new(x(1), x(2)) + beta*delta/(1.+(visits(x(1), x(2))/Na**.5)**.6) 
				
				theta_new(x(1), x(2)) = modulo(theta_new(x(1), x(2)) + &
						alpha/(1.+(visits(x(1), x(2))/Na**.5)**.3)*delta*elgb(a), 2.*pi)
	
				! save new occupation matrix
				visits(x_new(1), x_new(2)) = visits(x_new(1), x_new(2)) + 1
	
				occupation_new(x(1), x(2)) = occupation_new(x(1), x(2)) - 1
				occupation_new(x_new(1), x_new(2)) = occupation_new(x_new(1), x_new(2)) + 1
	
				! update state
				state_new(i,:) = x_new
	
			end do
			
			state = state_new
			v = v_new
			rewbar = rewbar_new
			occupation = occupation_new
	
			! after Trlx relaxation steps, sample the empirical density every 100 steps
			if ((ts > Trlx).and.(mod(ts-Trlx, 100).eq.1).and.(it == its)) then
				rho_av = rho_av + 1./float((T-Trlx)/100)*(occupation/(Na*dx**2.))
				v_av = v_av + 1./float((T-Trlx)/100)*v
			end if

			!	Print the average reward per walker at some time steps
			if (mod(ts,int(.01/dt))==0) then
				write(out_cost,*) ts*dt, avrew , rewbar
			end if
	
		end do
		write(out_cost,*) ''
		write(out_cost,*) ''


		! at the end of each iteration, save
		! ..states
		open(unit=in_states, file=trim(adjustl(in_name_states)), action="write", status="replace")
		! ..thetas
		open(unit=in_vt, file=trim(adjustl(in_name_vt)), action="write", status="replace")
		do i=1, Na
			write(in_states,*) state(i,:)
		end do
		write(in_vt,*) rewbar
		do i=-L,L
			do j=-L, L
				write(in_vt,*) v(i,j), theta(i,j)
			end do
		end do
		close(in_states)
		close(in_vt)
				
	end do
	! ---------------------- end of the learning -------------------------------


	!
	!	PRINT FINAL RESULTS
	!
	! --------------------------------------------------------------------------
	! print the average density and value (averaged over last T-Trlx --every 100-- steps)
	do i=-L, L
		do j=-L, L
			write(out_rho_av,*) i*dx, j*dx, rho_av(i,j)
		end do
		write(out_rho_av,*) ''
	end do
	radial_rho = axis_average_1s(rho_av)	! radial_average(rho_av)
	radial_val = axis_average_1s(v_av)		! radial_average(v_av)
	do i=0, L
		write(out_radial,*) i*dx, radial_rho(i), radial_val(i)
	end do
	
	do i=1, Ntraj
		close(out_traj(i))
	end do
	close(out_cost)
	close(out_theta)
	close(out_rho_value)
	close(out_rho_av)
	close(out_radial)

	write(100+proc,fmt="(a, i2, a, f7.4, a)") "Process ", proc, ", g = ", g, " --> finished!"

call MPI_Finalize ( ierr )

contains

	! cost function (function of points on the lattice)
	function q (x)
		implicit none
		integer, intent(in) :: x(2)
		real(dp) :: q
		
		q = .5*k*(x(1)**2. + x(2)**2.)*(Dx**2.)
	end function q

	! cost for collisions
	function collision (x)
		implicit none
		integer, intent(in) :: x(2)
		real(dp) :: collision
!		integer :: pairs
!		pairs = occupation(x(1),x(2))*(occupation(x(1),x(2)) - 1)
		collision = g*(Na-1)*dx**2.*(occupation(x(1),x(2))-1)
	end function collision


	subroutine choose_action (a)
		integer :: a

		call random_number(r)
			
		! choose action
		if (r .lt. p(0)) then
			a = 0 ! take action STAY
		else if (r .lt. sum(p(0:1))) then
			a = 1 ! take action UP
		else if (r .lt. sum(p(0:2))) then
			a = 2 ! take action RIGHT
		else if (r .lt. sum(p(0:3))) then
			a = 3 ! take action DOWN
		else
			a = 4 ! take action LEFT
		end if
	end subroutine choose_action


	function radial_average (f)
		real(dp) :: f(-L:L,-L:L)
		real(dp) :: radial_average(0:L)
		integer :: counter(0:L)
		integer :: i, j, ind
		counter = 0
		radial_average = 0
		do i=-L, L
			do j=-L, L
				ind = int(sqrt(i**2. + j**2.))
				if (ind <= L) then
					counter(ind) = counter(ind) + 1
					radial_average(ind) = radial_average(ind) + f(i,j)
				end if
			end do
		end do
		do i=0, L
			radial_average(i) = radial_average(i)/counter(i)
		end do
	end function radial_average


	function axis_average (f)
		real(dp) :: f(-L:L,-L:L)
		real(dp) :: axis_average(-L:L)
		integer :: i, j, ind
		axis_average = 0
		do i=-L, L
			axis_average(i) = (f(i,0)+f(0,i))/2.
		end do
	end function axis_average


	function axis_average_1s (f)
		real(dp) :: f(-L:L,-L:L)
		real(dp) :: axis_average_1s(0:L)
		integer :: i, j, ind
		axis_average_1s = 0
		do i=0, L
			axis_average_1s(i) = (f(i,0)+f(0,i)+f(-i,0)+f(0,-i))/4.
		end do
	end function axis_average_1s


	subroutine set_IO ()
		!
		!	OUTPUTS
		!
		! standard output
		out_form = "(a, i3, a3, f3.1, a1, i1,a1,i3)"
		write(out_name, out_form) "output/stdout_", 100+proc, "_N_", 1.*Na/(10**int(log10(1.*Na))) , "e", int(log10(1.*Na))
		print *, proc, "   ", trim(adjustl(out_name))
		open(100+proc, file=trim(adjustl(out_name)), action="write", status="replace")
	
		! create/clean output directories
		out_form = "(a, f4.2, a3, f3.1, a1, i1, a3, f6.4)"
		write(out_dir,out_form) "output/typlen_", d/U, &
								"/N_", 1.*Na/(10**int(log10(1.*Na))) , "e", int(log10(1.*Na)), &
								"/g_", g
	
		! backup old files
		call system("for i in $(ls "//trim(adjustl(out_dir))//"/*.dat); do mv -f $i $i'~'; rm -f $i; done")
		call system("for i in $(ls "//trim(adjustl(out_dir))//"/trajs/*.dat); do mv -f $i $i'~'; rm -f $i; done")
		
		! create directory and open files for the first Ntraj trajectories
		out_form = "(a, a5)"
		write(mk_dir,out_form) trim(adjustl(out_dir)) // "/trajs"
		call system("mkdir -p " // trim(adjustl(mk_dir)))
		out_form = "(a, i4, a4)"
		do i=1, Ntraj
			write(out_name,trim(out_form)) trim(adjustl(out_dir)) // "/trajs/traj_", i, ".dat"
			call replace_blanks(out_name,'0')
			open(unit=out_traj(i), file=trim(adjustl(out_name)), action="write", status="replace")
		end do
		
		! open file for the average cost
		out_form = "(a, f5.3, a4 )"
		write(out_name,out_form) trim(adjustl(out_dir)) // "/avg_cost_dx", dx, ".dat"
		call replace_blanks(out_name,'0')
		open(unit=out_cost, file=trim(adjustl(out_name)), action="write", status="replace")
		
		! open file for the drift field
		write(out_name,out_form) trim(adjustl(out_dir)) // "/drift_dx", dx,".dat"
		call replace_blanks(out_name,'0')
		open(unit=out_theta, file=trim(adjustl(out_name)), action="write", status="replace")
		
		! open file for the density and value on the lattice ..
		write(out_name,out_form) trim(adjustl(out_dir)) // "/rho_value_dx", dx, ".dat"
		call replace_blanks(out_name,'0')
		open(unit=out_rho_value, file=trim(adjustl(out_name)), action="write", status="replace")
		! .. and as radial average
		write(out_name,out_form) trim(adjustl(out_dir)) // "/radial_rho_value_dx", dx, ".dat"
		call replace_blanks(out_name,'0')
		open(unit=out_radial, file=trim(adjustl(out_name)), action="write", status="replace")
		
		! open file for the empirical density
		write(out_name,*) trim(adjustl(out_dir)) // "/rho_av.dat"
		open(unit=out_rho_av, file=trim(adjustl(out_name)), action="write", status="replace")
	
	
		!
		!	INPUT
		!
		! states and policy of last timestep of the previous simulation is the input
		! iostat gives 0 if the file exits, or something else if some error occurs
		!
		! file containing the positions of the walkers
		write(in_name_states,*) trim(adjustl(out_dir)) // "/previous_states.dat"
		open(unit=in_states, file=trim(adjustl(in_name_states)), action="read", status="old", iostat=states_error)
		! file containing the policy parameters
		write(in_name_vt,*) trim(adjustl(out_dir)) // "/previous_value_thetas.dat"
		open(unit=in_vt, file=trim(adjustl(in_name_vt)), action="read", status="old", iostat=vt_error)
	end subroutine set_IO

end program eikonal2d