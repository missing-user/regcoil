module init_surface_mod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine init_surface(lscreen_optin)

      use compute_offset_surface_mod
      use regcoil_variables, only: R0_plasma, nfp, volume_coil, &
                               ntheta_coil, nzeta_coil, nzetal_coil, &
                               theta_coil, zeta_coil, zetal_coil, &
                               r_coil, drdtheta_coil, drdzeta_coil, &
                               normal_coil, norm_normal_coil, area_coil, &
                               geometry_option_coil, R0_coil, a_coil,  &
                               separation, dtheta_coil, dzeta_coil, &
                               nescin_filename, mpol_coil, ntor_coil, &
                               nescin_filename, curpol

      use stel_kinds
      use stel_constants
      use omp_lib

      implicit none

      !character(*) :: nescin_filename
      !integer :: geometry_option, iflag, which_surface
      integer :: iflag, which_surface
      ! These next 2 arrays are not used now, but may be needed in the future for Merkel's regularization
      real(dp), dimension(:,:,:), allocatable :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2
      !real(dp), dimension(:,:), allocatable :: norm_normal, major_R_squared
      real(dp), dimension(:,:), allocatable :: major_R_squared
      real(dp) :: R0_to_use
      real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta
      real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
      real(dp) :: d2sinangledtheta2, d2cosangledtheta2, d2sinangle2dzeta2, d2cosangle2dzeta2
      integer :: i, itheta, izeta, im, in, mnmax_ws, istat = 0, iunit = 8
      real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_theta, delta_zeta, temp
      real(dp) :: r_tmp, rc_tmp, rs_tmp, zc_tmp, zs_tmp, factor

      ! variables to handle printing to the screen
      logical, optional :: lscreen_optin
      logical :: lscreen
    
      if(present(lscreen_optin)) then 
        lscreen = lscreen_optin
      else
        lscreen = .true.
      endif
    
      which_surface = 1
      
      ! Adding some checks to release previously allocated variables.
      ! This is because STELLOPT may call this function multiple times.
      ! if (allocated()) deallocate()
 

      if (allocated(theta_coil)) deallocate(theta_coil)
      allocate(theta_coil(ntheta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 1'

      if (allocated(zeta_coil)) deallocate(zeta_coil)
      allocate(zeta_coil(nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 2'

      if (allocated(zetal_coil)) deallocate(zetal_coil)
      allocate(zetal_coil(nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 3'

      do i = 1,ntheta_coil
         theta_coil(i) = twopi*(i-1.0_dp)/ntheta_coil
      end do

      do i = 1,nzeta_coil
         zeta_coil(i) = twopi/nfp*(i-1.0_dp)/nzeta_coil
      end do

      do i = 1,nzetal_coil
         zetal_coil(i) = twopi*(i-1.0_dp)/nzetal_coil
      end do

      dtheta_coil = theta_coil(2)-theta_coil(1)
      dzeta_coil  = zeta_coil(2)-zeta_coil(1)


      ! First dimension is the Cartesian component x, y, or z.
      if (allocated(r_coil)) deallocate(r_coil)
      allocate(r_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 4'

      if (allocated(drdtheta_coil)) deallocate(drdtheta_coil)
      allocate(drdtheta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 5'

      if (allocated(drdzeta_coil)) deallocate(drdzeta_coil)
      allocate(drdzeta_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 6'

      if (allocated(normal_coil)) deallocate(normal_coil)
      allocate(normal_coil(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 7'

      r_coil = 0
      drdtheta_coil = 0
      drdzeta_coil = 0

      ! We do not use the 2nd derivatives presently, but these placeholders are here in case we need them in the future:
      if (allocated(d2rdtheta2)) deallocate(d2rdtheta2)
      allocate(d2rdtheta2(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 8'
      if (allocated(d2rdthetadzeta)) deallocate(d2rdthetadzeta)
      allocate(d2rdthetadzeta(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 9'
      if (allocated(d2rdzeta2)) deallocate(d2rdzeta2)
      allocate(d2rdzeta2(3,ntheta_coil,nzetal_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 10'

      d2rdtheta2 = 0
      d2rdthetadzeta = 0
      d2rdzeta2 = 0

      if (geometry_option_coil==3 .or. geometry_option_coil == 4) then
         if (lscreen) print *,"  Reading coil surface from nescin file ",trim(nescin_filename)

         call read_nescin(nescin_filename, r_coil, drdtheta_coil, drdzeta_coil, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, &
              ntheta_coil, nzetal_coil, theta_coil, zetal_coil, .false.)
      end if


      select case (geometry_option_coil)
      case (0,1)
         ! Torus with circular cross-section

         if (lscreen) print *,"  Building a plain circular torus."

         if (geometry_option_coil==0) then
            R0_to_use = R0_plasma
         else
            R0_to_use = R0_coil
         end if

         do itheta = 1,ntheta_coil
            angle = theta_coil(itheta)
            sinangle = sin(angle)
            cosangle = cos(angle)
            dsinangledtheta = cosangle
            dcosangledtheta = -sinangle
            d2sinangledtheta2 = -sinangle
            d2cosangledtheta2 = -cosangle
            do izeta = 1,nzetal_coil
               angle2 = zetal_coil(izeta)
               sinangle2 = sin(angle2)
               cosangle2 = cos(angle2)
               dsinangle2dzeta = cosangle2
               dcosangle2dzeta = -sinangle2
               d2sinangle2dzeta2 = -sinangle2
               d2cosangle2dzeta2 = -cosangle2

               r_coil(1,itheta,izeta) = (R0_to_use + a_coil * cosangle) * cosangle2
               r_coil(2,itheta,izeta) = (R0_to_use + a_coil * cosangle) * sinangle2
               r_coil(3,itheta,izeta) = a_coil * sinangle

               drdtheta_coil(1,itheta,izeta) = (a_coil * dcosangledtheta) * cosangle2
               drdtheta_coil(2,itheta,izeta) = (a_coil * dcosangledtheta) * sinangle2
               drdtheta_coil(3,itheta,izeta) = a_coil * dsinangledtheta

               drdzeta_coil(1,itheta,izeta) = (R0_to_use + a_coil * cosangle) * dcosangle2dzeta
               drdzeta_coil(2,itheta,izeta) = (R0_to_use + a_coil * cosangle) * dsinangle2dzeta
               !drdzeta(3,itheta,izeta) = 0, so no equation needed for it here.

! This next bit is remarked out since we don't need 2nd derivatives presently
!!$               if (transfer_matrix_option==2 .and. which_surface == 1) then
!!$                  d2rdu2(1,itheta,izeta) = a * d2cosangledu2 * cosangle2
!!$                  d2rdu2(2,itheta,izeta) = a * d2cosangledu2 * sinangle2
!!$                  d2rdu2(3,itheta,izeta) = a * d2sinangledu2
!!$
!!$                  d2rdudzeta(1,itheta,izeta) = a * dcosangledu * dcosangle2dv
!!$                  d2rdudzeta(2,itheta,izeta) = a * dcosangledu * dsinangle2dv
!!$                  !d2rdudzeta(3,itheta,izeta) = 0, so no equation needed for it here.
!!$
!!$                  d2rdv2(1,itheta,izeta) = (R0_to_use + a * cosangle) * d2cosangle2dv2
!!$                  d2rdv2(2,itheta,izeta) = (R0_to_use + a * cosangle) * d2sinangle2dv2
!!$                  !d2rdv2(3,itheta,izeta) = 0, so no equation needed for it here.
!!$               end if
            end do
         end do

      case (2,4)

!!$         if (transfer_matrix_option==2) then
!!$            stop "Error! This geometry_option_middle is not yet implemented for transfer_matrix_option=2."
!!$         end if

         if (geometry_option_coil==2) then
            if (lscreen) print "(a,f10.4,a)","   Constructing a surface offset from the plasma by ",separation," meters."
            ! print "(a,f10.4,a)","   Constructing a surface offset from the plasma by ",separation," meters."
         else
            if (lscreen) print "(a,f10.4,a)","   Constructing a surface offset from the nescin surface by ",separation," meters."
         end if

         ! Finite differences to use:
         ! (Numerical Recipes suggests (machine epsilon)^(1/3)
         delta_theta = 1e-5;
         delta_zeta = 1e-5;
         ! Trick from Numerical Recipes for improved accuracy:
         temp = 1.0_dp + delta_theta
         delta_theta = temp - 1.0_dp
         temp = 1.0_dp + delta_zeta
         delta_zeta = temp - 1.0_dp
 

         !$OMP PARALLEL

         !$OMP MASTER
         if (lscreen) print *,"  Number of OpenMP threads:",omp_get_num_threads()
         !$OMP END MASTER

         !$OMP DO PRIVATE(x_new,y_new,z_new,x_old,y_old,z_old)
         do itheta = 1,ntheta_coil
            do izeta = 1,nzetal_coil

               ! Compute r:
               call compute_offset_surface_xyz_of_thetazeta(theta_coil(itheta),zetal_coil(izeta),x_new,y_new,z_new,separation)
               r_coil(1,itheta,izeta) = x_new
               r_coil(2,itheta,izeta) = y_new
               r_coil(3,itheta,izeta) = z_new

               ! Compute dr/du:
               call compute_offset_surface_xyz_of_thetazeta(theta_coil(itheta)-delta_theta,zetal_coil(izeta),x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_thetazeta(theta_coil(itheta)+delta_theta,zetal_coil(izeta),x_new,y_new,z_new,separation)
               drdtheta_coil(1,itheta,izeta) = (x_new-x_old)/(2*delta_theta)
               drdtheta_coil(2,itheta,izeta) = (y_new-y_old)/(2*delta_theta)
               drdtheta_coil(3,itheta,izeta) = (z_new-z_old)/(2*delta_theta)

               ! Compute dr/dv:
               call compute_offset_surface_xyz_of_thetazeta(theta_coil(itheta),zetal_coil(izeta)-delta_zeta,x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_thetazeta(theta_coil(itheta),zetal_coil(izeta)+delta_zeta,x_new,y_new,z_new,separation)
               drdzeta_coil(1,itheta,izeta) = (x_new-x_old)/(2*delta_zeta)
               drdzeta_coil(2,itheta,izeta) = (y_new-y_old)/(2*delta_zeta)
               drdzeta_coil(3,itheta,izeta) = (z_new-z_old)/(2*delta_zeta)

            end do

         end do
         !$OMP END DO
         !$OMP END PARALLEL                  
         call write_nescin  ! write the offset winding surface
         
      case (3)

         !call write_nescin  ! write the offset winding surface
         ! Nothing to do - we already read the nescin file.

      case default
         print *,"Invalid setting for geometry_option_coil: ",geometry_option_coil
         stop
      end select

      ! Evaluate cross product:
      normal_coil(1,:,:) = drdzeta_coil(2,:,:) * drdtheta_coil(3,:,:) - drdtheta_coil(2,:,:) * drdzeta_coil(3,:,:)
      normal_coil(2,:,:) = drdzeta_coil(3,:,:) * drdtheta_coil(1,:,:) - drdtheta_coil(3,:,:) * drdzeta_coil(1,:,:)
      normal_coil(3,:,:) = drdzeta_coil(1,:,:) * drdtheta_coil(2,:,:) - drdtheta_coil(1,:,:) * drdzeta_coil(2,:,:)

      if (allocated(norm_normal_coil)) deallocate(norm_normal_coil)
      allocate(norm_normal_coil(ntheta_coil, nzeta_coil),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error! init_surface_mod 11'
      norm_normal_coil = sqrt(normal_coil(1,:,1:nzeta_coil)**2 + normal_coil(2,:,1:nzeta_coil)**2 &
                               +  normal_coil(3,:,1:nzeta_coil)**2)

      area_coil = nfp * dtheta_coil * dzeta_coil * sum(norm_normal_coil)

      ! Compute coil surface volume using \int (1/2) R^2 dZ dzeta.
      ! These quantities will be evaluated on the half theta grid, which is the natural grid for dZ,
      ! but we will need to interpolate R^2 from the full to half grid.
      allocate(major_R_squared(ntheta_coil,nzetal_coil))
      major_R_squared = r_coil(1,:,:)*r_coil(1,:,:) + r_coil(2,:,:)*r_coil(2,:,:)
      ! First handle the interior of the theta grid:
      volume_coil = sum((major_R_squared(1:ntheta_coil-1,:) + major_R_squared(2:ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r_coil(3,2:ntheta_coil,:)-r_coil(3,1:ntheta_coil-1,:))) ! dZ
      ! Add the contribution from the ends of the theta grid:
      volume_coil = volume_coil + sum((major_R_squared(1,:) + major_R_squared(ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
           * (r_coil(3,1,:)-r_coil(3,ntheta_coil,:))) ! dZ
      volume_coil = abs(volume_coil * dzeta_coil / 2) ! r includes all nfp periods already, so no factor of nfp needed.
      deallocate(major_R_squared)
      if (lscreen) print "(a,es10.3,a,es10.3,a)"," Coil surface area:",area_coil," m^2. Volume:",volume_coil," m^3."

    end subroutine init_surface

end module init_surface_mod
