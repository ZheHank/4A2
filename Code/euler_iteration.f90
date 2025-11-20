
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(:,:), allocatable :: mass_i, flux_i
      real, dimension(:,:), allocatable :: mass_j, flux_j
      integer :: i, j, ni, nj

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERTED
      allocate(mass_i(ni,nj-1), mass_j(ni-1,nj))
      allocate(flux_i(ni,nj-1), flux_j(ni-1,nj))
      mass_i = 0.5 * (g%rovx(1:ni,1:nj-1) + g%rovx(1:ni,2:nj)) * g%lx_i(1:ni,1:nj-1)&
            + 0.5 * (g%rovy(1:ni,1:nj-1) + g%rovy(1:ni,2:nj)) * g%ly_i(1:ni,1:nj-1)
      mass_j = 0.5 * (g%rovx(1:ni-1,1:nj) + g%rovx(2:ni,1:nj)) * g%lx_j(1:ni-1,1:nj) &
            + 0.5 * (g%rovy(1:ni-1,1:nj) + g%rovy(2:ni,1:nj)) * g%ly_j(1:ni-1,1:nj)
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

      ! modify to use starting values for Runge-Kutta -Added
      g%ro = g%ro_start
!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERTED
      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro,g%dro)

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
      flux_i = mass_i * 0.5 * (g%hstag(1:ni,1:nj-1) + g%hstag(1:ni,2:nj))
      flux_j = mass_j * 0.5 * (g%hstag(1:ni-1,1:nj) + g%hstag(2:ni,1:nj))

! modify to use starting values for Runge-Kutta -Added
      g%roe = g%roe_start
!     Update the internal energy with enthalpy fluxes
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g%area,g%roe,g%droe)

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
      flux_i = mass_i * 0.5 * (g%vx(1:ni,1:nj-1) + g%vx(1:ni,2:nj)) &
            + 0.5 * (g%p(1:ni,1:nj-1) + g%p(1:ni,2:nj)) * g%lx_i(1:ni,1:nj-1)
      flux_j = mass_j * 0.5 * (g%vx(1:ni-1,1:nj) + g%vx(2:ni,1:nj)) &
            + 0.5 * (g%p(1:ni-1,1:nj) + g%p(2:ni,1:nj)) * g%lx_j(1:ni-1,1:nj)
!     Update the x-momentum with momentum flux

! modify to use starting values for Runge-Kutta -Added
      g%rovx = g%rovx_start
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx,g%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
      flux_i = mass_i * 0.5 * (g%vy(1:ni,1:nj-1) + g%vy(1:ni,2:nj)) &
            + 0.5 * (g%p(1:ni,1:nj-1) + g%p(1:ni,2:nj)) * g%ly_i(1:ni,1:nj-1)
      flux_j = mass_j * 0.5 * (g%vy(1:ni-1,1:nj) + g%vy(2:ni,1:nj)) & 
            + 0.5 * (g%p(1:ni-1,1:nj) + g%p(2:ni,1:nj)) * g%ly_j(1:ni-1,1:nj)

! modify to use starting values for Runge-Kutta -Added
      g%rovy = g%rovy_start
!     Update the y-momentum with momentum flux
!     INSERT
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy,g%drovy)

!     Add artificial viscosity by smoothing all of the primary flow variables

      call smooth_array(av,g%ro,g%corr_ro,av%ro_ref)
      call smooth_array(av,g%roe,g%corr_roe,av%roe_ref)
      call smooth_array(av,g%rovx,g%corr_rovx,av%rov_ref)
      call smooth_array(av,g%rovy,g%corr_rovy,av%rov_ref)
      

      end subroutine euler_iteration


