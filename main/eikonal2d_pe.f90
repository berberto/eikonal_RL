!
!    Module with constants and functions
!
module consts_funcs

    use kinds, only: dp
    use constants
    use mpi

    implicit none

    real(dp), parameter :: k = 1.                ! curvature harmonic potential
    real(dp), parameter :: D = .0001             ! diffusion constant
    real(dp), parameter :: U = .0002             ! speed
    real(dp), parameter :: Dx = .05              ! lattice spacing
    real(dp), parameter :: Dt = 1.               ! time discretization
    real(dp), parameter :: R = 10.               ! square domain (-R,R)x(-R,R)
    integer, parameter :: L = ceiling(R/dx)      ! lattice (-L,L)x(-L,L)
    integer, parameter :: T    = int(60000./dt)! time steps per iteration

    real(dp), parameter :: eps = .5*U*Dt/Dx
    real(dp), parameter :: sigma = D*Dt/(Dx**2.)

    private :: sigma, eps, R

    public :: Dx, Dt, k, L, T, U, D

contains

    subroutine check_params ()
        if((sigma < eps).or.(sigma > .5)) then
            print *, "Check parameters, Dx and Dt"
            call exit(1)
        end if
    end subroutine check_params


    ! policy function
    function policy (theta)
        real(dp) :: theta
        real(dp) :: policy(0:4)

        policy(0) = 1.- 4.*sigma           ! STAY
        policy(1) = sigma + eps*sin(theta) ! UP
        policy(2) = sigma + eps*cos(theta) ! RIGHT
        policy(3) = sigma - eps*sin(theta) ! DOWN
        policy(4) = sigma - eps*cos(theta) ! LEFT
    end function policy

end module consts_funcs



!
!    MAIN PROGRAM
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
    integer :: Na                        ! number of agents

    ! random number
    integer, allocatable :: seed(:)
    real(dp) :: r

    ! learning variables
    real(dp) :: p(0:4)                        ! policy
    integer :: i, j, a, ts                    ! running indices
    integer :: x(2), x_new(2)                 ! temporary arrays for single-agent states
    integer, allocatable :: state(:,:)        ! first index is the agent, second is the lattice position

    integer :: start
    integer :: move(0:4, 2)         ! first index is the action, second is lattice displacement
    real(dp) :: theta(-L:L, -L:L)   ! table of the parameter of the policy
    real(dp) :: rew, avrew, ret

    ! output
    integer :: out_val
    character(len=64) :: out_form, out_name, out_dir

    ! check that policy is well defined
    call check_params

    !
    !    Parameters from command line
    !
    if (iargc() < 1) then
        print *, "Enter number of particles"
        call exit (1)
    end if

    call getarg (1,arg)
    read(arg,*) Na
    
    !
    !    Initialize MPI environment
    !
call MPI_Init ( ierr )
    call MPI_Comm_size(MPI_COMM_WORLD, Nproc, error)
    call MPI_Comm_rank(MPI_COMM_WORLD, proc, error)

    if((Nproc > 20).and.(proc .eq. 0)) then
        print *, "WARNING: Running on more than 20 cores..."
!        print *, "Run with maximum 20 processors:"
!        print *, "  mpirun -np <#proc> ./eikonal <#agents (parallel runs)>"
!        call exit(1)
    end if

    out_val = 100+proc
    start = 2*proc

    allocate(state(Na,2))

    !
    !    Set the INPUT and OUTPUT files
    !
    ! create/clean output directories
    out_form = "(a, f4.2, a10, f6.4)"
    write(out_dir,out_form) "output/typlen_", d/U, "/p_ev_free"
    call system("mkdir -p " // trim(adjustl(out_dir)))

    ! open output files: one per value of g and per starting point
    out_form = "(a, f4.2, a4)"
    write(out_name, out_form) "/diff_value_", start*Dx, ".dat"
    open(unit=out_val, file=trim(adjustl(out_dir))//trim(adjustl(out_name)), action="write", status="replace")

    print *, proc, trim(adjustl(out_dir))//trim(adjustl(out_name))
!    call MPI_Finalize(ierr)
!    stop

    ! inizialize random number generator
    call seed_from_urandom(seed)

    ! define the actions (displacements on the lattice)
    move(0,:) = (/ 0, 0 /)    ! STAY in place
    move(1,:) = (/ 0, 1 /)    ! move UP
    move(2,:) = (/ 1, 0 /)    ! move RIGHT
    move(3,:) = (/ 0, -1 /)   ! move DOWN
    move(4,:) = (/ -1, 0 /)   ! move LEFT


    !
    !    INITIALIZATION OF STATES AND POLICY
    !
    ! --------------------------------------------------------------------------
    ! initialization of the position of the Na independend walkers
    do i=1, Na
        state(i,:) = (/ start, 0 /)
    end do
    
    ! central drift
    do i=-L, L
        do j=-L, L
            if (j > 0) then
                theta(i,j) = pi + acos((i-.5)/sqrt((i-.5)**2. + (j-.5)**2.))
            else
                theta(i,j) = pi - acos((i-.5)/sqrt((i-.5)**2. + (j-.5)**2.))
            end if
        end do
    end do
   

    !
    !    START THE POLICY EVALUATION
    !
    ! --------------------------------------------------------------------------
    ! Take T steps for each iteration, for its iterations
    ret = 0.
    do ts=0, T

        !
        !   For each particle take action, calculate reward
        !   and average reward over all the particles
        !
        avrew = 0.
        do i=1, Na
            x = state(i,:)
            p = policy(theta(x(1), x(2)))

            ! choose single-agent action according to the policy
            call choose_action(a)

            ! take action: move in the opposite direction if going outside the domain
            x_new(:) = x(:) + move(a,:)
            if ((abs(x_new(1)) .gt. L) .or. (abs(x_new(2)) .gt. L)) then
                x_new(:) = x(:) - move(a,:)
            end if

            ! save new state
            state(i,:) = x_new

            ! calculate average reward per particle
            avrew = avrew - q(x_new)/Na

        end do

        ! update return
        ret = ret + avrew

        ! Print the average reward per walker at some time steps
        if (mod(ts,int(20/dt))==0) then
            write(out_val,"(3es14.4)") ts*dt, ret
        end if

    end do
    ! ------------------ end of the policy-evaluation --------------------------

    close(out_val)
call MPI_Finalize ( ierr )

contains

    ! cost function (function of points on the lattice)
    function q (x)
        implicit none
        integer, intent(in) :: x(2)
        real(dp) :: q

        q = .5*k*((x(1)-.5)**2. + (x(2)-.5)**2.)*(Dx**2.)
    end function q


    subroutine choose_action (a)
        integer :: a

        call random_number(r)

        ! choose action
        if (r .lt. p(0)) then
            a = 0 ! STAY
        else if (r .lt. sum(p(0:1))) then
            a = 1 ! UP
        else if (r .lt. sum(p(0:2))) then
            a = 2 ! RIGHT
        else if (r .lt. sum(p(0:3))) then
            a = 3 ! DOWN
        else
            a = 4 ! LEFT
        end if

    end subroutine choose_action




end program eikonal2d
