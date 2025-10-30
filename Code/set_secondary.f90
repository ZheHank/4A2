      
      subroutine set_secondary(av,g)

!     This subroutine calculates the secondary flow variables from the primary
!     ones at every node

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g

!     Define any further variables you may need
!     INSERT
      real, dimension(:,:), allocatable :: rho, roe, rovx, rovy
      real, dimension(:,:), allocatable :: p, hstag
      integer :: ni, nj

!     The primary flow variables are "ro", "roe", "rovx" and "rovy", these are 
!     the conserved quantities within the Euler equations. Write the code to
!     calculate the secondary flow variables which are the velocity components
!     "vx" and "vy", the static pressure "p" and the stagnation enthalpy
!     "hstag". These are needed at every timestep, there is no need for any 
!     loops as the operations can be performed elementwise, although you may
!     wish to define some intermediate variables to improve readability.
!     INSERTED
      ni = g%ni
      nj = g%nj
      allocate(rho(ni,nj), roe(ni,nj), rovx(ni,nj), rovy(ni,nj))
      allocate(p(ni,nj), hstag(ni,nj))
      rho = g%ro
      roe = g%roe
      rovx = g%rovx
      rovy = g%rovy

      g%vx = rovx / rho
      g%vy = rovy / rho
      p = (roe - 0.5 * rho * (g%vx**2 + g%vy**2)) * (av%gam - 1.0)
      g%p = p
      hstag = (roe + p) / rho
      g%hstag = hstag

      end subroutine set_secondary


