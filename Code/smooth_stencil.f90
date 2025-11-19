
      module smooth_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop,corr)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_2, prop_avg_4
      integer :: ni, nj
      real, parameter :: fcorr = 0.9, flinear = 0.25
      real, dimension(size(prop,1),size(prop,2)), intent(inout) :: corr
      real, dimension(size(prop,1),size(prop,2)) :: corr_total

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
!     INSERT
      prop_avg_2(2:ni-1,2:nj-1) = 0.25 * (prop(1:ni-2,2:nj-1) + prop(3:ni,2:nj-1) + &
                           prop(2:ni-1,1:nj-2) + prop(2:ni-1,3:nj))
      

!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
!     INSERT
      prop_avg_2(1,2:nj-1) = (prop(1,1:nj-2) + prop(1,3:nj) + 2.0 * prop(2,2:nj-1) - prop(3,2:nj-1)) / 3.0
      prop_avg_2(ni,2:nj-1) = (prop(ni,1:nj-2) + prop(ni,3:nj) + 2.0 * prop(ni-1,2:nj-1) - prop(ni-2,2:nj-1)) / 3.0
      prop_avg_2(2:ni-1,1) = (prop(1:ni-2,1) + prop(3:ni,1) + 2.0 * prop(2:ni-1,2) - prop(2:ni-1,3)) / 3.0
      prop_avg_2(2:ni-1,nj) = (prop(1:ni-2,nj) + prop(3:ni,nj) + 2.0 * prop(2:ni-1,nj-1) - prop(2:ni-1,nj-2)) / 3.0

!     The corner values are not currently smoothed
      prop_avg_2([1,ni],[1,nj]) = prop([1,ni],[1,nj])
! Now repeat the process to get a four-point average.(Lagrange interpolation)
      prop_avg_4(:,:) = prop_avg_2(:,:)
      prop_avg_4(3:ni-2, 3:nj-2) = (-1.0/6.0) * prop(1:ni-4, 3:nj-2) + (2.0/3.0) * prop(2:ni-3, 3:nj-2) + &
                              (2.0/3.0) * prop(4:ni-1, 3:nj-2) + (-1.0/6.0) * prop(5:ni, 3:nj-2)

      prop_avg_4(3:ni-2, 3:nj-2) = prop_avg_4(3:ni-2, 3:nj-2) + (-1.0/6.0) * prop(3:ni-2, 1:nj-4) + &
                              (2.0/3.0) * prop(3:ni-2, 2:nj-3) + (2.0/3.0) * prop(3:ni-2, 4:nj-1) + &
                              (-1.0/6.0) * prop(3:ni-2, 5:nj)

      ! Average of both directions
      prop_avg_4(3:ni-2, 3:nj-2) = prop_avg_4(3:ni-2, 3:nj-2) / 2.0

      ! Edge values for the four-point average (for top, bottom, left, right edges)
      ! Bottom edge (i = 3 to ni-2, j = 1)
      prop_avg_4(3:ni-2, 1) = (-1.0/6.0) * prop(1:ni-4, 1) + (2.0/3.0) * prop(2:ni-3, 1) + &
                              (2.0/3.0) * prop(4:ni-1, 1) + (-1.0/6.0) * prop(5:ni, 1)

      ! Top edge (i = 3 to ni-2, j = nj)
      prop_avg_4(3:ni-2, nj) = (-1.0/6.0) * prop(1:ni-4, nj) + (2.0/3.0) * prop(2:ni-3, nj) + &
                              (2.0/3.0) * prop(4:ni-1, nj) + (-1.0/6.0) * prop(5:ni, nj)

      ! Left edge (i = 1, j = 3 to nj-2)
      prop_avg_4(1, 3:nj-2) = (-1.0/6.0) * prop(1, 1:nj-4) + (2.0/3.0) * prop(1, 2:nj-3) + &
                              (2.0/3.0) * prop(1, 4:nj-1) + (-1.0/6.0) * prop(1, 5:nj)
      ! Right edge (i = ni, j = 3 to nj-2)
      prop_avg_4(ni, 3:nj-2) = (-1.0/6.0) * prop(ni, 1:nj-4) + (2.0/3.0) * prop(ni, 2:nj-3) + &
                              (2.0/3.0) * prop(ni, 4:nj-1) + (-1.0/6.0) * prop(ni, 5:nj)

!     Corners for the four-point average can be set to the original values or handled
!     with a simple average of available neighbors.
      prop_avg_4(1,1) = prop(1,1)
      prop_avg_4(1,nj) = prop(1,nj)
      prop_avg_4(ni,1) = prop(ni,1)
      prop_avg_4(ni,nj) = prop(ni,nj)
!     Combine the two averages to get a smoother result


!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
!     INSERT
!      prop = (1.0 - av%sfac) * prop + av%sfac * prop_avg

!     Now with defferred correction for spatial accuracy
      corr_total = fcorr * (prop - flinear * prop_avg_2 - (1.0 - flinear) * prop_avg_4)

      corr = 0.99 * corr + 0.01 * corr_total
      prop = (1.0 - av%sfac) * prop + av%sfac * (flinear * prop_avg_2 + (1.0 - flinear) * prop_avg_4 + corr)

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil


