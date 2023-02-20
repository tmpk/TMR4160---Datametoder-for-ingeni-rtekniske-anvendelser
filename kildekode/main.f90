program main
	! Program for solving the 2-D Poisson equation											
	!     d^2            d^2
	!    ----- u(x,y) + ----- u(x,y) = g(x,y)
	!	  dx^2			 dx^2
	! on the square unit grid, by discretizing the equation and using the iterative Jacobi- or Gauss-Seidel
	! methods. User may choose grid step size (identical in x- and y-direction) and 
	! boundary conditions through console I/O. Also allows a preset two-choice of g(x,y)
	! for the Neumann b.c. Results are written to file.
	!
	! For the Dirichlet b.c, g(x,y) = 1, while for Neumann b.c either
	! g(x,y) = 12 - 12x - 12y', or
	! g(x,y) = (6 - 12x) * (3y^2 - 2y^3) + (3x^2 - 2x^3)*(6-12y)'
	! 
	! In the case of Jacobi iteration, the equation is solved by explicitly applying the five-point
	! formula.
	!
	! If using Gauss-Seidel iteration, the equations is solved similarly to the Jacobi procedure, OR 
	! a system of equations Au=b is formulated and solved.
	! 
	! 
	! Date/version: 30.04.2019/1.0
	
	implicit none
	double precision, parameter :: tolerance = 0.000001			! error tolerance
	integer :: i,j												! loop variables
	integer :: iter												! number of iterations 
	character(len=3) :: solver									! specifies solver algorithm
	integer numx, numy											! number of nodes is numx*numy		
	double precision, dimension(:,:), allocatable :: u			! matrix containing solution u
	double precision, dimension(:,:), allocatable :: g			! matrix containing rhs of Poisson eq.
	logical :: dirichlet, neumann_g1, neumann_g2				! boundary conditions
	logical :: explicit											! specifies variant of GS method
	real :: h													! grid step size
	character :: keypress										! user input 
	double precision time_start, time_end						! helper variables for measuring time
	
	write(*,*)'This program solves the 2D Poisson equation'
	write(*,*)
	write(*,*) '    d^2            d^2'
	write(*,*) '   ----- u(x,y) + ----- u(x,y) = g(x,y)'
	write(*,*) '    dx^2           dx^2'
	write(*,*)
	write(*,*)'for the unknown function u on the square unit grid.'
	write(*,*)
	
	! prompt user for choosing tolerance, grid step size, boundary conditions and g(x,y) (if choosing Neumann bc)
100	call user_input(h, dirichlet,neumann_g1,neumann_g2,solver,explicit)
	numx = 1.0/h + 1
	numy=numx
	
	allocate(u(numx,numy))
	
	allocate(g(numx,numy))
	call build_g(g,dirichlet,neumann_g1,neumann_g2, numx, numy, h)
	
	write(*,*)'Working...'
	
	call cpu_time(time_start)
	if (solver=='jm') then
		call poisson_solver_jm(u,g,numx,numy,tolerance,h,iter,dirichlet)
	else
		call poisson_solver_gsm(u,g,numx,numy,tolerance,iter,dirichlet,h,explicit)
	end if	
	call cpu_time(time_end)
	
	write(*,*)'System was solved in:',(time_end - time_start),'seconds.'   

!-------------print u to console if it is not too large
	if (numx<=11) then
		write(*,*)'Computed solution u is'
		write(*,*)
		do i=1,numx
    		write(*,"(100g15.5)") ( u(i,j), j=1,numx )
		end do
		write(*,*)
	end if
!--------------------


!-------------print number of iterations 
	write(*,*) 'Number of iterations:',(iter)
!--------------------
	
! create and open new .dat-file for storing solution
! the solution is stored in format
!-------------
! numx numy
! u(1,1)    
! u(2,1)	
! ...
! u(numx,1)    
! u(1,2)	
! u(2,2)    
! ...
! ...
! u(numx, numy)  
!---------------
! Write solution to "solution.dat". If file does not exist it is created. 
! If it does exist, it is overwritten.
	open(100, file='solution.dat', status='replace')
	write(100,*)numx, numy
	do j=1,numy
		do i=1,numx
			write(100,*)u(i,j)
		end do
	end do
	close(100)
	write(*,*)'Results are written to file "solution.dat".'
	
	! deallocate allocated variables	
	deallocate(g)
	deallocate(u)

	write(*,*)
	write(*,'(A)', ADVANCE='NO')' Re-run program? (NOTE: this overwrites the "solution.dat"-file) y/n: ' 
	read(*,*) keypress
	if (keypress=='y') then
		write(*,*)
		go to 100
	end if

	write(*,*) 'Program terminating.'
	write(*,*)
	
end program

function r8mat_rms ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RMS returns the root mean square of data stored as an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the data whose RMS is desired.
!
!    Output, real ( kind = 8 ) R8MAT_RMS, the root mean square of A.
!
  implicit none

  integer m
  integer n

  double precision a(m,n)
  double precision r8mat_rms

  r8mat_rms = sqrt ( sum ( a(1:m,1:n)**2 ) / real ( m * n) )

  return
  
end 

subroutine user_input(h, dirichlet, neumann_g1, neumann_g2,solver,explicit)
! This routine prompts the user for choosing step size, boundary conditions and solver algorithm,
! and outputs variables corresponding to these
!
! 
! Date/version: 30.04.2019/1.0

	implicit none
	real, intent(OUT) :: h										! step size
	logical, intent(OUT) :: dirichlet, neumann_g1, neumann_g2	! boundary conditions
	logical, intent(OUT) :: explicit							! specifies GS method variant
	character(len=3), intent(OUT) :: solver						! solver algorithm
	character :: keypress1, keypress2, keypress3, keypress4		! user input
												! 							

101 write(*,'(A)',ADVANCE='NO')' Set the step size h (h must be in (0,1) and give equally spaced grid nodes): '
	read(*,'(f10.9)') h
	if (h>=1.0 .or. h<=0.0 .or. mod(1.0/h,1.0)/=0) then
		write(*,*)'Invalid input given...'
		go to 101
	else
		write(*,*)'Step size set to:',(h)
		write(*,*)
	end if
	
	write(*,*)'Boundary condition choices are either'
	write(*,*) ' 1) Dirichlet b.c (u is known on the boundary)'
	write(*,*) ' 2) Neumann b.c (du/dn is known on the boundary)'
	 do while (.not. (keypress1 == '1' .or. keypress1 == '2'))	
		write(*,'(A)',ADVANCE='NO')' Which boundary condition should be applied? 1/2: '
		read(*,*) keypress1
	end do
	
	
	if (keypress1 == '1') then
		dirichlet = .true.
		write(*,*)'Dirichlet boundary condition is applied.'
	else 
		dirichlet = .false.
		write(*,*)'Neumann boundary condition is applied.'
		write(*,*)
		write(*,*)'Possible choices for function g(x,y) are'
		write(*,*) '1) g(x,y) = 12 - 12x - 12y'
		write(*,*) '2) g(x,y) = (6 - 12x) * (3y^2 - 2y^3) + (3x^2 - 2x^3)*(6-12y)'
		do while (.not. (keypress2 == '1' .or. keypress2 == '2'))
				write(*,'(A)',ADVANCE='NO')' Which function g(x,y) should be used? 1/2: '	
				read(*,*) keypress2
		end do
		
		if (keypress2=='1') then
			neumann_g1 = .true.
			neumann_g2 = .false.
			write(*,*)'g(x,y) set to: g(x,y) = 12 - 12x - 12y.'
		else	
			neumann_g2 = .true.
			neumann_g1 = .false.
			write(*,*)'g(x,y) set to: g(x,y) = (6 - 12x) * (3y^2 - 2y^3) + (3x^2 - 2x^3)*(6-12y).'
		end if
		
	end if
	
	write(*,*)
	
	write(*,*)'Choose solver algorithm:'
	write(*,*) '1) Jacobi method'
	write(*,*) '2) Gauss-Seidel method'
	do while (.not. (keypress3 == '1' .or. keypress3 == '2'))
		write(*,'(A)',ADVANCE='NO')' Which solver should be used? 1/2: '	
		read(*,*) keypress3	
	end do
	
	if (keypress3=='1') then
		solver = 'jm'
		write(*,*)'Jacobi method chosen as solver algorithm.'
		write(*,*)
	else	
		solver = 'gsm'
		write(*,*)'Gauss-Seidel method chosen as solver algorithm.'
		write(*,*)
		write(*,*)'Use' 
		write(*,*)' 1) general formulation of GS algorithm, or'
		write(*,*)' 2) special formulation.'
		do while (.not. (keypress4 == '1' .or. keypress4 == '2'))
		write(*,'(A)',ADVANCE='NO')' Which algorithm variant should be used? 1/2: '	
		read(*,*) keypress4	
		end do
		if (keypress4 == '2') then
			explicit = .true.
		else
			explicit = .false.
		end if		
	end if

end subroutine

subroutine build_g(g, dirichlet, neumann_g1, neumann_g2, numx, numy, h)
! This routine creates a 2D matrix g that holds the values of the right hand side of the
! Poisson equation such that g(x,y)=g(i*h,j*h); where x,y are coordinates of the nodes,
! i,j are the indices of g and h is the grid step size.
! 
! In addition, if the Dirichlet boundary condition is applied, the routine
! sets the points on the matrix's boundary to the known value of u there.
! 
! 
! Date/version: 30.04.2019/1.0

	implicit none
	logical, intent(IN) :: dirichlet, neumann_g1, neumann_g2	! boundary conditions
	integer, intent(IN) :: numx, numy							! numx*numy is number of nodes
	real, intent(IN) :: h										! step size
	double precision, dimension(numx,numy), intent(OUT) :: g	! matrix g
	double precision :: x, y									! x- and y- coords
	integer i,j													! iteration variables
	
	if (dirichlet .eqv. .true.) then 
		g = 1
		! set the known boundary point values:
		do i=1,numx
			g(i,1)=0.25*((h*(i-1))**2)						
			g(i,numy)=0.25*((h*(i-1))**2+(h*(numy-1))**2)	
			g(1,i)=g(i,1)									
			g(numx,i)=g(i,numy)
		end do
	else 
		do j=1, numy
			x=(j-1)*h
			do i=1,numx
				y=(i-1)*h
				if (neumann_g1 .eqv. .true.) then
					g(i,j) = 12-12*x-12*y
				else if (neumann_g2 .eqv. .true.) then
					g(i,j) = (6-12*x)*(3*y**2-2*y**3) + (3*x**2-2*x**3)*(6-12*y)
				end if	
			end do
		end do
	end if
	
end subroutine

subroutine poisson_solver_jm(u,g,numx,numy,tolerance,h,iter,dirichlet)
! This routine solves the 2D Poisson equation by making the call
! to the appropriate subroutine, according to which boundary condition
! is active.
!
! 
! Date/version: 01.05.2019/1.0

	implicit none
	integer, intent(IN) :: numx,numy							! number of nodes is numx*numy
	real, intent(IN) :: h										! grid step size
	double precision, dimension(numx,numy), intent(IN) :: g		! matrix containing rhs of Poisson eq.
	double precision, intent(IN) :: tolerance					! error tolerance
	double precision, dimension(numx,numy), intent(OUT) :: u	! matrix containing solution u
	integer, intent(OUT) :: iter								! number of iterations 
	logical, intent(IN) :: dirichlet

	if (dirichlet .eqv. .true.) then
		call solve_dirichlet(u,g,numx,numy,tolerance,h,iter)
	else
		call solve_neumann(u,g,numx,numy,tolerance,h,iter)
	end if

end subroutine


subroutine solve_dirichlet(u, g, numx, numy, tolerance, h,iter)
! Solves the 2D Poisson equation with Dirichlet boundary condition by
! using the Jacobi iteration. Solution is returned in matrix u, and number of 
! iterations in variable iter
!
!
! Date/version: 30.04.2019/1.0

	implicit none
	integer, intent(IN) :: numx, numy							! number of nodes is numx*numy
	real, intent(IN) :: h										! step size
	double precision, dimension(numx,numy), intent(IN) :: g		! rhs of Poisson eq
	double precision, intent(IN) :: tolerance					! error tolerance
	double precision, dimension(numx,numy), intent(OUT) :: u	! solution matrix
	integer, intent(OUT) :: iter								! number of iterations
	integer :: i, j												! loop variables
	double precision, dimension(numx,numy) :: udiff				! differences between current and prev. iteration
	double precision, dimension(numx,numy) :: u_new				! vector with new iterate
	double precision :: diff									! difference between rms of current and prev. iteration
	logical :: done												! controls if were done iterating
	double precision r8mat_rms									! function for calculating rms								

	u = 0			! initial guess for solution

	iter = 1					
	done = .false.
	do while (.not. done)
		do j=1, numy
			do i=1, numx
			if (i==1 .or. j==1 .or. i==numx .or. j==numy) then
				! on boundary, so value is known 
				u_new(i,j)=g(i,j)								
			else
				! use five-point formula for unknown values on interior
				u_new(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-g(i,j)*h**2) 
			end if
			end do
		end do
		udiff = u_new-u
		diff = r8mat_rms(numx,numy,udiff)
		if (diff <= tolerance) then 
			done = .true.
		else
			u=u_new
		end if
		iter=iter+1
	end do

end subroutine

subroutine solve_neumann(u,g,numx,numy,tolerance,h,iter)
! Solves the Poisson equation on the square unit grid,
! with Neumann boundary condition du/dn=0, using the Jacobi method.
! Returns the solution in the matrix u and number of iterations performed in iter.
! 
! 
! Date/version: 30.04.2019/1.0

	implicit none
	integer, intent(IN) :: numx, numy							! number of nodes is numx*numy
	double precision, dimension(numx,numy), intent(IN) :: g 	! rhs of Poisson eq
	double precision, intent(IN) :: tolerance					! error tolerance
	real, intent(IN) :: h										! step size
	integer, intent(OUT) :: iter								! number of iterations
	double precision, dimension(numx,numy), intent(OUT) :: u	! solution matrix
	double precision :: diff									! difference between rms of current and prev. iteration
	double precision, dimension(numx,numy) :: udiff		 		! differences between current and prev. iteration
	double precision, dimension(numx,numy) :: u_new				! matrix with most recent iterates
	integer :: i, j												! loop variables
	logical :: done												! controls whether we're done iterating
	double precision r8mat_rms									! function for returning rms of data

	u = 0		! initial guess for solution

	done = .false.
	iter = 1
	do while (.not. done)
		do j=1,numy
			do i=1,numx
			! formula for u(i,j) is different for  
			! 1) the 4 corners
			! 2) the 4 boundaries excluding corner points
			! 3) interior
			! => 9 different formulas for u(i,j):
				if (j==1) then
					if (i==1) then								! lower left corner
						u_new(i,j)=0
					else if (i==numx) then						! upper left corner
						u_new(i,j)=0.25*(2*u(i-1,j)+2*u(i,j+1)-g(i,j)*h**2)
					else										! left boundary excluding corners
						u_new(i,j)=0.25*(u(i+1,j)+u(i-1,j)+2*u(i,j+1)-g(i,j)*h**2)
					end if
				else if (j==numy) then		
					if (i==1) then								! lower right corner
						u_new(i,j)=0.25*(2*u(i+1,j)+2*u(i,j-1)-g(i,j)*h**2)
					else if (i==numx) then						! upper right corner
						u_new(i,j)=0.25*(2*u(i-1,j)+2*u(i,j-1)-g(i,j)*h**2)
					else										! right boundary excluding corners
						u_new(i,j)=0.25*(u(i+1,j)+u(i-1,j)+2*u(i,j-1)-g(i,j)*h**2)
					end if
				else if (i==1 .and. j>1 .and. j<numy) then		! lower boundary excluding corners
					u_new(i,j)=0.25*(2*u(i+1,j)+u(i,j-1)+u(i,j+1)-g(i,j)*h**2)
				else if (i==numx .and. j>1 .and. j<numy) then	! upper boundary excluding corners
					u_new(i,j)=0.25*(2*u(i-1,j)+u(i,j-1)+u(i,j+1)-g(i,j)*h**2)
				else											! interior
					u_new(i,j)=0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-g(i,j)*h**2)
				end if
			end do
		end do
		udiff = u_new-u
		diff = r8mat_rms(numx,numy,udiff)
		if (diff <= tolerance) then 
			done = .true.
		else
			u=u_new
		end if
		iter = iter+1
	end do

end subroutine

subroutine poisson_solver_gsm( u,g,numx,numy,tolerance,iter,dirichlet,h,explicit )
! This routine solves the 2D Poisson equation using the Gauss-Seidel method. It
! makes the appropriate call to a subroutine depending on which variant of the GS method
! should be used. The algorithm used is either one using the explicit formulas for u(i,j) or
! a more general matrix formulation. In the latter case, this routine must also 
! set up the system Au=b by making calls to subroutine.
!
! 
! Date/version 01.05.2019/1.0

	implicit none
	integer, intent(IN) :: numx, numy							! number of nodes is numx*numy
	logical, intent(IN) :: explicit								! specifies variant of GS method
	double precision, dimension(numx,numy), intent(OUT) :: u	! matrix containing solution u
	double precision, dimension(numx,numy), intent(IN) :: g		! matrix containing rhs of Poisson eq.
	double precision, dimension(:,:), allocatable :: A			! matrix A in Au=b
	double precision, dimension(:), allocatable :: b			! vector b in Au=b
	integer :: n												! number of unknown variables 
	double precision, intent(IN) :: tolerance
	logical, intent(IN) :: dirichlet
	integer, intent(OUT) :: iter
	real, intent(IN) :: h
	
	if (explicit .eqv. .true.) then
		call gsm_solver_explicit(u, g, numx, numy, tolerance, h,iter,dirichlet)
	else
		if (dirichlet .eqv. .true.) then
			n = numx - 2
		else
			n = numx
		end if
		
		allocate(A(n**2,n**2))
		call build_A(A,n,dirichlet)
	
		allocate(b(n**2))
		call build_b(b,g,n,h,dirichlet,numx)
	
		call gsm_solver(u,A,b,g,numx,numy,tolerance,iter,n,dirichlet)
	
		! deallocate allocated variables
		deallocate(A)
		deallocate(b)
	end if
	
end subroutine

subroutine gsm_solver_explicit(u, g, numx, numy, tolerance, h,iter, dirichlet)
! This routine solves the 2D Poisson equation by the iterative Gauss-Seidel method,
! and returns the solution in the matrix u and the number of iterations performed in the variable iter.
! It uses the explicit formulas for u(i,j), and essentially copies the code in solve_dirichlet and solve_neumann,
! replacing u with u_new at the appropriate places in each formula.
!
! 
! Date/version: 02.05.2019/1.0


	implicit none
	integer, intent(IN) :: numx, numy							! number of nodes is numx*numy
	real, intent(IN) :: h										! step size
	double precision, dimension(numx,numy), intent(IN) :: g		! rhs of Poisson eq
	double precision, intent(IN) :: tolerance					! error tolerance
	logical, intent(IN) :: dirichlet							! boundary condition
	double precision, dimension(numx,numy), intent(OUT) :: u	! solution matrix
	integer, intent(OUT) :: iter								! number of iterations
	integer :: i, j												! loop variables
	double precision, dimension(numx,numy) :: udiff				! differences between current and prev. iteration
	double precision, dimension(numx,numy) :: u_new				! vector with new iterate
	double precision :: diff									! difference between rms of current and prev. iteration
	logical :: done												! controls if were done iterating
	double precision r8mat_rms									! function for calculating rms								

	u = 0			! initial guess for solution
	iter = 1					
	done = .false.
	
	if (dirichlet .eqv. .true.) then
		do while (.not. done)
			do j=1, numy
				do i=1, numx
				if (i==1 .or. j==1 .or. i==numx .or. j==numy) then
					! on boundary, so value is known 
					u_new(i,j)=g(i,j)								
				else
					! use five-point formula for unknown values on interior
					u_new(i,j) = 0.25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1)+u(i,j+1)-g(i,j)*h**2) 
				end if
				end do
			end do
			udiff = u_new-u
			diff = r8mat_rms(numx,numy,udiff)
			if (diff <= tolerance) then 
				done = .true.
			else
				u=u_new
			end if
			iter=iter+1
		end do
	else
		do while (.not. done)
			do j=1,numy
				do i=1,numx
				! formula for u(i,j) is different for  
				! 1) the 4 corners
				! 2) the 4 boundaries excluding corner points
				! 3) interior
				! => 9 different formulas for u(i,j):
					if (j==1) then
						if (i==1) then								! lower left corner
							u_new(i,j)=0
						else if (i==numx) then						! upper left corner
							u_new(i,j)=0.25*(2*u_new(i-1,j)+2*u(i,j+1)-g(i,j)*h**2)
						else										! left boundary excluding corners
							u_new(i,j)=0.25*(u(i+1,j)+u_new(i-1,j)+2*u(i,j+1)-g(i,j)*h**2)
						end if
					else if (j==numy) then		
						if (i==1) then								! lower right corner
							u_new(i,j)=0.25*(2*u(i+1,j)+2*u_new(i,j-1)-g(i,j)*h**2)
						else if (i==numx) then						! upper right corner
							u_new(i,j)=0.25*(2*u_new(i-1,j)+2*u_new(i,j-1)-g(i,j)*h**2)
						else										! right boundary excluding corners
							u_new(i,j)=0.25*(u(i+1,j)+u_new(i-1,j)+2*u_new(i,j-1)-g(i,j)*h**2)
						end if
					else if (i==1 .and. j>1 .and. j<numy) then		! lower boundary excluding corners
						u_new(i,j)=0.25*(2*u(i+1,j)+u_new(i,j-1)+u(i,j+1)-g(i,j)*h**2)
					else if (i==numx .and. j>1 .and. j<numy) then	! upper boundary excluding corners
						u_new(i,j)=0.25*(2*u_new(i-1,j)+u_new(i,j-1)+u(i,j+1)-g(i,j)*h**2)
					else											! interior
						u_new(i,j)=0.25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1)+u(i,j+1)-g(i,j)*h**2)
					end if
				end do
			end do
			udiff = u_new-u
			diff = r8mat_rms(numx,numy,udiff)
			if (diff <= tolerance) then 
				done = .true.
			else
				u=u_new
			end if
			iter = iter+1
		end do
	end if
	
end subroutine


subroutine build_A(A, n, dirichlet)
! This routine builds the matrix A in Au=b. The A matrix is formed from matrices D and the identity matrix.
! The size and structure of A and D depends on whether the Dirichlet or Neumann b.c is active.
!
! 
! Date/version: 30.04.2019/1.0
	
	integer, intent(in) :: n									! number of unknowns in system
	logical, intent(in) :: dirichlet							! Dirichlet boundary condition
	double precision, dimension(n**2,n**2), intent(out) :: A	! system matrix A
	double precision, dimension(n,n) :: eye, D					! helper matrices for constructing A

	! build identity matrix eye
	do j=1,n
		do i=1,n
			if (i==j) then
				eye(i,j)=1
			else
				eye(i,j)=0
			end if
		end do
	end do
	
	! build D matrix
	do j=1,n
		do i=1,n
			if (i==j) then								
				D(i,j)=4								
				if (i==1) then							
					if (dirichlet .eqv. .false.) then
						D(i,j+1) = -2										
					else 								
						D(i,j+1) = -1					
					end if
					D(i+1,j) = -1						
				else if (i==n) then 
					if (dirichlet .eqv. .false.) then 
						D(i,j-1)=-2
					else 
						D(i,j-1)=-1
					end if
					D(i-1,j)=-1
				else
					D(i+1,j)=-1
					D(i,j+1)=-1
				end if
			end if
		end do
	end do	
	
	! build A by inserting D and eye at the appropriate
	! indices. 
	do i=0,n*(n-1),n
		A(i+1:i+n,i+1:i+n)=D
		if (i==0) then
			if (dirichlet .eqv. .false.) then
				A(i+1:i+n,i+n+1:i+2*n)=-2*eye
			else
				A(i+1:i+n,i+n+1:i+2*n)=-eye
			end if
		else if (i==n*(n-1)) then
			if (dirichlet .eqv. .false.) then
				A(i+1:i+n,i-n+1:i)=-2*eye
			else
				A(i+1:i+n,i-n+1:i)=-eye
			end if
		else
			A(i+1:i+n,i-n+1:i)=-eye
			A(i+1:i+n,i+n+1:i+2*n)=-eye
		end if
	end do
	
	! apply Neumann b.c to A if it is active
	if (dirichlet .eqv. .false.) then
		do i=2,n*n
			A(1,i)=0
		end do
	end if
		
end subroutine

subroutine build_b(b,g,n,h,dirichlet,numx)
! This routine creates the b vector in Au=b. If the Dirichlet b.c is applied, the known boundary values of u
! must be included in b. The other values in b are simply retrieved from g 
!
! 
! Date/version: 30.04.2019/1.0

	implicit none
	integer, intent(in) :: n								! number of unknowns in each row
	integer, intent(in) :: numx								! grid size is numx*numx
	double precision, dimension(numx,numx), intent(in) :: g	! rhs of Poisson eq. and known boundary points if Dirichlet b.c
	real, intent(in) :: h									! step size
	logical, intent(in) :: dirichlet						! Dirichlet b.c
	double precision, dimension(n**2), intent(out) :: b		! vector b in Au=b
	integer :: counter, k, l, i, j							! helper variables for counting, looping/indexing	
	
	counter=1
	
	if (dirichlet .eqv. .true.) then
	! must remember to include known values at the boundaries in b.
	! These values are stored in g.
		k=2
		l=numx-1
		! loop over interior points as these are the unknowns:
		do j=k,l
			do i=k,l
				if (i==k) then
					if (j==k) then
						b(counter)=g(i-1,j)+g(i,j-1)-g(i,j)*h**2
					else if (j==l) then
						b(counter)=g(i-1,j)+g(i,j+1)-g(i,j)*h**2
					else
						b(counter)=g(i-1,j)-g(i,j)*h**2
					end if
				else if (i==l) then
					if (j==k) then
						b(counter)=g(i+1,j)+g(i,j-1)-g(i,j)*h**2
					else if (j==l) then
						b(counter)=g(i+1,j)+g(i,j+1)-g(i,j)*h**2
					else
						b(counter)=g(i+1,j)-g(i,j)*h**2
					end if
				else if (j==k .and. i>k .and. i<l) then
					b(counter)=g(i,j-1)-g(i,j)*h**2
				else if (j==l .and. i>k .and. i<l) then
					b(counter)=g(i,j+1)-g(i,j)*h**2
				else
					b(counter)=-g(i,j)*h**2
				end if
				counter=counter+1
			end do
		end do
	else
	! values in b are found from g
		k=1
		l=numx
	! loop over all points on grid:
		do j=k,l
			do i=k,l
				b(counter)=-g(i,j)*h**2
				counter=counter+1
			end do
		end do
	end if
	
	! apply Neumann b.c to b if it is active
	if (dirichlet .eqv. .false.) then
		b(1)=0
	end if
	
end subroutine

subroutine gsm_solver(u,A,b,g,numx,numy, tolerance, iter,n, dirichlet)
! This routine solves the 2D Poisson equation by the iterative Gauss-Seidel method,
! and returns the solution in the matrix u and the number of iterations performed in the variable iter.
! It uses the more general formulation of the algorithm:
! https://www.cfd-online.com/Wiki/Gauss-Seidel_method
!
! 
! Date/version: 30.04.2019/1.0

	implicit none
	integer, intent(in) :: numx, numy							! number of nodes is numx*numy
	integer, intent(in) :: n									! number of unknowns
	double precision, dimension (n**2,n**2), intent(in) :: A	! matrix A in Au=b
	double precision, dimension(numx,numy), intent(in) :: g  	! matrix g with rhs of Poisson equation
	double precision, dimension(n**2), intent(in) :: b			! vector b in Au=b
	double precision, intent(in) :: tolerance					! tolerance used for terminating iterations
	logical, intent(in) :: dirichlet							! dirichlet boundary condition
	double precision, dimension(numx,numy), intent(out) :: u	! solution matrix
	integer, intent(out) :: iter								! number of iterations
	double precision :: left_sum, right_sum		 ! left_sum/right_sum is sum of A(i,j)*x(j) to the LEFT/RIGHT  of diagonal of A				
	logical :: done												! logical for checking convergence
	double precision :: diff									! difference between rms of new and old iterate
	double precision, dimension(n*n) :: x, x_new				! vectors containing old and new iterate
	double precision, dimension(n*n) :: x_diff					! vector containing difference between iterates
	double precision r8mat_rms									! function for returning rms of data
	integer :: i,j,counter										! iteration and counter variables 

	x = 0
	x_new = 0
	done = .false.
	iter = 1
	
	do while (.not. done)
		do i = 1,n*n
			left_sum = 0
			right_sum = 0
			do j = 1,i-1
				left_sum = left_sum+A(i,j)*x_new(j)
			end do
			do j = i+1,n*n
				right_sum = right_sum+A(i,j)*x(j)
			end do
			x_new(i) = (1.0/A(i,i))*(b(i)-left_sum-right_sum)
		end do
		x_diff = x_new-x
		diff = r8mat_rms(n*n,1,x_diff)
		if (diff <= tolerance) then 
			done = .true.
		else
			x=x_new
		end if
		iter = iter+1
	end do

	! output the solution as a matrix with correct boundary values 
	counter = 1
	do j=1,numx
		do i=1,numx
			if (dirichlet .eqv. .true.) then
				if (i==1 .or. i==numx .or. j==1 .or. j==numx) then
					u(i,j)=g(i,j)
				else 
					u(i,j) = x_new(counter)
					counter = counter+1
				end if
			else
				u(i,j) = x_new(counter)
				counter = counter+1
			end if
		end do
	end do

end subroutine
