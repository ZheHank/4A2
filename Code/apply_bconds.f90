      
      subroutine apply_bconds(av,g,bcs)

!     This subroutine applies both the inlet and outlet boundary conditions, as
!     it modifies both the primary and secondary flow variables they must be
!     calculated first

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs

!     Declare the other variables you need here
!     INSERT
      real, dimension(:), allocatable :: rho_in, t_in, v_in, e_in, hstag_in
      real, dimension(:), allocatable :: p_in, roe_in, vx_in, vy_in, rovx_in, rovy_in
      real :: alpha_rad

!     At the inlet boundary the change in density is driven towards "rostag",
!     which is then used to obtain the other flow properties to match the
!     specified stagnation pressure, temperature and flow angle. 

!     To help prevent instabilities forming at the inlet boundary condition the 
!     changes in inlet density are relaxed by a factor "rfin" normally set to 
!     0.25 but it can be reduced further.

!     It is also worth checking if "ro" is greater than "rostag" and limiting 
!     the values to be slightly less than "rostag". This can prevent the solver 
!     crashing during severe transients.
      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)
      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

!     Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!     "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!     "hstag(1,:)"
!     INSERTED
      allocate(rho_in(size(bcs%ro)), t_in(size(bcs%ro)), v_in(size(bcs%ro)), e_in(size(bcs%ro)), hstag_in(size(bcs%ro)))
      allocate(p_in(size(bcs%ro)), roe_in(size(bcs%ro)), vx_in(size(bcs%ro)), vy_in(size(bcs%ro)), rovx_in(size(bcs%ro)), rovy_in(size(bcs%ro)))
      bcs%rostag = bcs%pstag / (av%rgas * bcs%tstag)
      rho_in = bcs%ro
      t_in = bcs%tstag * (rho_in / bcs%rostag)**(av%gam - 1.0)
      v_in = sqrt(2.0 * av%cp * (bcs%tstag - t_in))
      e_in = av%cv * t_in + 0.5 * v_in**2
      hstag_in = av%cp * bcs%tstag
      p_in = rho_in * av%rgas * t_in
      roe_in = rho_in * e_in
      alpha_rad = bcs%alpha * (3.14159265358979323846 / 180.0)
      vx_in = v_in * cos(alpha_rad)
      vy_in = v_in * sin(alpha_rad)
      rovx_in = rho_in * vx_in
      rovy_in = rho_in * vy_in

      g%p(1,:) = p_in
      g%rovx(1,:) = rovx_in
      g%rovy(1,:) = rovy_in
      g%roe(1,:) = roe_in
      g%vx(1,:) = vx_in
      g%vy(1,:) = vy_in
      g%hstag(1,:) = hstag_in

!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
!     INSERTED
      g%p(av%ni,:) = bcs%p_out

      end subroutine apply_bconds


