#:include 'macros.fpp'

module m_additional_forcing

    use m_derived_types
    use m_global_parameters
    use m_ibm
    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_additional_forcing_module, s_finalize_additional_forcing_module, s_compute_periodic_forcing, &
        & s_update_controllers

    ! forcing params
    type(scalar_field), allocatable, dimension(:) :: q_periodic_force
    real(wp)                                      :: avg_coeff
    real(wp)                                      :: spatial_rho, spatial_rhou, spatial_rhoe
    real(wp), allocatable, dimension(:)           :: rho_window, rhou_window, rhoe_window
    real(wp)                                      :: sum_rho, sum_rhou, sum_rhoe
    real(wp)                                      :: phase_rho, phase_rhou, phase_rhoe
    integer                                       :: window_fill

    $:GPU_DECLARE(create='[q_periodic_force, avg_coeff]')
    $:GPU_DECLARE(create='[spatial_rho, spatial_rhou, spatial_rhoe, phase_rho, phase_rhou, phase_rhoe]')

    ! control params
    real(wp)                            :: rho_avg_loc, rhou_avg_loc, cs_avg_loc
    real(wp), allocatable, dimension(:) :: err_u_hist
    real(wp), allocatable, dimension(:) :: rho_wdw_cntrl, u_wdw_cntrl, cs_wdw_cntrl, Vp_wdw_cntrl
    real(wp)                            :: rho_sum_cntrl, u_sum_cntrl, cs_sum_cntrl, Vp_sum_cntrl
    integer                             :: wdw_fill_cntrl

    $:GPU_DECLARE(create='[rho_avg_loc, rhou_avg_loc, cs_avg_loc]')

contains

    subroutine s_initialize_additional_forcing_module

        integer  :: i
        real(wp) :: domain_vol

        ! total cartesian domain volume
        domain_vol = (x_domain%end - x_domain%beg)*(y_domain%end - y_domain%beg)*(z_domain%end - z_domain%beg)

        ! coefficient used for phase averages
        avg_coeff = 1._wp/(domain_vol*fluid_volume_fraction)
        $:GPU_UPDATE(device='[avg_coeff]')

        if (periodic_forcing) then
            @:ALLOCATE(q_periodic_force(1:num_dims+2))
            do i = 1, num_dims + 2
                @:ALLOCATE(q_periodic_force(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(q_periodic_force(i))
            end do

            ! initialization of parameters
            window_fill = 0

            @:ALLOCATE(rho_window(forcing_window))
            @:ALLOCATE(rhou_window(forcing_window))
            @:ALLOCATE(rhoe_window(forcing_window))

            rho_window = 0.0_wp
            rhou_window = 0.0_wp
            rhoe_window = 0.0_wp

            sum_rho = 0.0_wp
            sum_rhou = 0.0_wp
            sum_rhoe = 0.0_wp

            phase_rho = 0._wp
            phase_rhou = 0._wp
            phase_rhoe = 0._wp

            if (forcing_wrt .and. proc_rank == 0) then
                open (unit=102, file='forcing.bin', status='replace', form='unformatted', access='stream', action='write')
            end if
        end if

        if (particle_control) then
            @:ALLOCATE(err_u_hist(4))
            err_u_hist = 0._wp

            @:ALLOCATE(rho_wdw_cntrl(cntrl_p%window_size))
            @:ALLOCATE(u_wdw_cntrl(cntrl_p%window_size))
            @:ALLOCATE(cs_wdw_cntrl(cntrl_p%window_size))
            @:ALLOCATE(Vp_wdw_cntrl(cntrl_p%window_size))
            rho_wdw_cntrl = 0._wp
            u_wdw_cntrl = 0._wp
            cs_wdw_cntrl = 0._wp
            Vp_wdw_cntrl = 0._wp

            rho_sum_cntrl = 0._wp
            u_sum_cntrl = 0._wp
            cs_sum_cntrl = 0._wp
            Vp_sum_cntrl = 0._wp

            wdw_fill_cntrl = 0

            if (forcing_wrt .and. proc_rank == 0) then
                open (unit=103, file='controls.bin', status='replace', form='unformatted', access='stream', action='write')
            end if
        end if

    end subroutine s_initialize_additional_forcing_module

    !> compute the space and time average of quantities, compute the periodic forcing terms described in Khalloufi and Capecelatro
    subroutine s_compute_periodic_forcing(rhs_vf, q_cons_vf, q_prim_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        integer, intent(in)                                    :: t_step
        real(wp)                                               :: spatial_rho_glb, spatial_rhou_glb, spatial_rhoe_glb
        real(wp)                                               :: dVol, rho, f_rhou
        integer                                                :: window_loc
        integer                                                :: i, j, k, l

        ! zero spatial averages
        spatial_rho = 0._wp
        spatial_rhou = 0._wp
        spatial_rhoe = 0._wp

        $:GPU_UPDATE(device='[spatial_rho, spatial_rhou, spatial_rhoe]')

        ! compute spatial averages
        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[spatial_rho, spatial_rhou, spatial_rhoe]]', reductionOp='[+]', private='[l, &
                            & rho, dVol]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then
                        rho = 0._wp
                        do l = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k)
                        end do
                        dVol = dx(i)*dy(j)*dz(k)
                        spatial_rho = spatial_rho + (rho*dVol)  ! rho
                        spatial_rhou = spatial_rhou + (q_cons_vf(eqn_idx%cont%end + mom_f_idx)%sf(i, j, k)*dVol)  ! rho*u
                        spatial_rhoe = spatial_rhoe + ((q_cons_vf(eqn_idx%E)%sf(i, j, &
                                                       & k) - 0.5_wp*rho*(q_prim_vf(eqn_idx%mom%beg)%sf(i, j, &
                                                       & k)**2 + q_prim_vf(eqn_idx%mom%beg + 1)%sf(i, j, &
                                                       & k)**2 + q_prim_vf(eqn_idx%mom%beg + 2)%sf(i, j, k)**2))*dVol)  ! rho*e
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[spatial_rho, spatial_rhou, spatial_rhoe]')

        ! reduction sum across entire domain
        call s_mpi_allreduce_sum(spatial_rho, spatial_rho_glb)
        call s_mpi_allreduce_sum(spatial_rhou, spatial_rhou_glb)
        call s_mpi_allreduce_sum(spatial_rhoe, spatial_rhoe_glb)

        spatial_rho_glb = spatial_rho_glb*avg_coeff
        spatial_rhou_glb = spatial_rhou_glb*avg_coeff
        spatial_rhoe_glb = spatial_rhoe_glb*avg_coeff

        ! update time average window location
        window_loc = 1 + mod(t_step, forcing_window)

        ! update time average sum
        sum_rho = sum_rho - rho_window(window_loc) + spatial_rho_glb
        sum_rhou = sum_rhou - rhou_window(window_loc) + spatial_rhou_glb
        sum_rhoe = sum_rhoe - rhoe_window(window_loc) + spatial_rhoe_glb

        ! update window arrays
        rho_window(window_loc) = spatial_rho_glb
        rhou_window(window_loc) = spatial_rhou_glb
        rhoe_window(window_loc) = spatial_rhoe_glb

        ! update number of time samples
        if (window_fill < forcing_window) window_fill = window_fill + 1

        ! compute phase averages
        phase_rho = sum_rho/real(window_fill, wp)
        phase_rhou = sum_rhou/real(window_fill, wp)
        phase_rhoe = sum_rhoe/real(window_fill, wp)
        $:GPU_UPDATE(device='[phase_rho, phase_rhou, phase_rhoe]')

        ! compute periodic forcing terms for mass, momentum, energy
        $:GPU_PARALLEL_LOOP(collapse=3, private='[l, rho, f_rhou]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then
                        rho = 0._wp
                        do l = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k)
                        end do
                        ! continuity forcing
                        q_periodic_force(1)%sf(i, j, k) = (rho_inf_ref - phase_rho)*forcing_dt

                        ! momentum forcing
                        f_rhou = (rho_inf_ref*u_inf_ref - phase_rhou)*forcing_dt
                        do l = 1, num_dims
                            q_periodic_force(1 + l)%sf(i, j, k) = 0._wp !q_prim_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, k)*q_periodic_force(1)%sf(i, j, k)
                        end do

                        q_periodic_force(1 + mom_f_idx)%sf(i, j, k) = q_periodic_force(1 + mom_f_idx)%sf(i, j, k) + f_rhou

                        ! energy forcing
                        q_periodic_force(2 + num_dims)%sf(i, j, k) = (P_inf_ref*gammas(1) - phase_rhoe)*forcing_dt & 
                            & + q_prim_vf(eqn_idx%mom%beg + mom_f_idx - 1)%sf(i, j, k)*f_rhou 
                            ! & + q_cons_vf(eqn_idx%E)%sf(i, j, k)*q_periodic_force(1)%sf(i, j, k)/rho
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! add the forcing terms to the RHS
        $:GPU_PARALLEL_LOOP(collapse=3, private='[l, rho]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then
                        rho = 0._wp
                        do l = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k)
                        end do
                        do l = 1, num_fluids
                            rhs_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k) = rhs_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, &
                                   & k) + q_cons_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k)*q_periodic_force(1)%sf(i, j, k)/rho  ! continuity
                        end do
                        ! energy
                        rhs_vf(eqn_idx%E)%sf(i, j, k) = rhs_vf(eqn_idx%E)%sf(i, j, k) + q_periodic_force(2 + num_dims)%sf(i, j, k)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (t_step > forcing_start) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[l]')
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == 0) then
                            do l = 1, num_dims
                                rhs_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, k) = rhs_vf(eqn_idx%mom%beg + l - 1)%sf(i, j, &
                                       & k) + q_periodic_force(1 + l)%sf(i, j, k)  ! momentum
                            end do
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (forcing_wrt .and. proc_rank == 0) then
            ! print *, 'FORCING:', spatial_rho_glb, spatial_rhou_glb, spatial_rhoe_glb
            write (102) spatial_rho_glb, spatial_rhou_glb, spatial_rhoe_glb
            flush (102)
        end if

    end subroutine s_compute_periodic_forcing

    subroutine s_update_controllers(t_step, q_cons_vf, q_prim_vf)

        integer, intent(in) :: t_step
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp) :: mu, Dp, gamma 
        real(wp) :: rho, dVol, pres
        real(wp) :: rho_avg, rhou_avg, u_avg, cs_avg
        real(wp) :: u_star_rel, Mach, u_rel
        real(wp) :: err_u, err_M, d_err_u
        real(wp) :: Vp_avg_loc, Vp_avg
        integer  :: ib_local
        integer :: window_loc
        integer :: i, j, k, l

        ! these need to be moved to case file/setup
        mu = 1._wp/fluid_pp(1)%Re(1)
        Dp = 0.1_wp ! hardcoded currently 
        gamma = 1._wp/fluid_pp(1)%gamma + 1._wp

        ! intialize
        Vp_avg_loc = 0._wp
        rho_avg_loc = 0._wp
        rhou_avg_loc = 0._wp
        cs_avg_loc = 0._wp

        do i = 1, num_local_ibs
            ib_local = local_ib_patch_ids(i)
            Vp_avg_loc = Vp_avg_loc + patch_ib(ib_local)%vel(mom_f_idx)
        end do

        ! Global sum over the unique owner partition.
        call s_mpi_allreduce_sum(Vp_avg_loc, Vp_avg)

        Vp_avg = Vp_avg / real(num_gbl_ibs, wp)

        ! get averages of density, momentum, soundspeed
        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[rho_avg_loc, rhou_avg_loc, cs_avg_loc]]', reductionOp='[+]', private='[l, rho, dVol, pres]', copyin='[gamma]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then
                        rho = 0._wp
                        do l = 1, num_fluids
                            rho = rho + q_cons_vf(eqn_idx%cont%beg + l - 1)%sf(i, j, k)
                        end do
                        dVol = dx(i)*dy(j)*dz(k)
                        rho_avg_loc = rho_avg_loc + (rho*dVol)
                        rhou_avg_loc = rhou_avg_loc + (q_cons_vf(eqn_idx%mom%beg + mom_f_idx - 1)%sf(i, j, k)*dVol)
                        pres = q_prim_vf(eqn_idx%E)%sf(i, j, k)
                        cs_avg_loc = cs_avg_loc + (sqrt(gamma*pres/rho)*dVol)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_UPDATE(host='[rho_avg_loc, rhou_avg_loc, cs_avg_loc]')

        ! reduction sum across entire domain
        call s_mpi_allreduce_sum(rho_avg_loc, rho_avg)
        call s_mpi_allreduce_sum(rhou_avg_loc, rhou_avg)
        call s_mpi_allreduce_sum(cs_avg_loc, cs_avg)

        rho_avg = rho_avg*avg_coeff
        rhou_avg = rhou_avg*avg_coeff
        cs_avg = cs_avg*avg_coeff
        u_avg = rhou_avg/rho_avg

        ! time average over window
        window_loc = 1 + mod(t_step, cntrl_p%window_size)

        rho_sum_cntrl = rho_sum_cntrl - rho_wdw_cntrl(window_loc) + rho_avg
        u_sum_cntrl = u_sum_cntrl - u_wdw_cntrl(window_loc) + u_avg
        cs_sum_cntrl = cs_sum_cntrl - cs_wdw_cntrl(window_loc) + cs_avg
        Vp_sum_cntrl = Vp_sum_cntrl - Vp_wdw_cntrl(window_loc) + Vp_avg

        rho_wdw_cntrl(window_loc) = rho_avg
        u_wdw_cntrl(window_loc) = u_avg
        cs_wdw_cntrl(window_loc) = cs_avg
        Vp_wdw_cntrl(window_loc) = Vp_avg

        if (wdw_fill_cntrl < cntrl_p%window_size) wdw_fill_cntrl = wdw_fill_cntrl + 1

        rho_avg = rho_sum_cntrl/real(wdw_fill_cntrl, wp)
        u_avg = u_sum_cntrl/real(wdw_fill_cntrl, wp)
        cs_avg = cs_sum_cntrl/real(wdw_fill_cntrl, wp)
        Vp_avg = Vp_sum_cntrl/real(wdw_fill_cntrl, wp)
        ! done time averaging

        u_rel = u_avg - Vp_avg

        u_star_rel = cntrl_p%Re_tgt*mu/(rho_avg*Dp)
        Mach = u_rel/cs_avg

        err_u = u_star_rel - u_rel
        err_M = cntrl_p%M_tgt - Mach

        err_u_hist(1) = err_u_hist(2)
        err_u_hist(2) = err_u_hist(3)
        err_u_hist(3) = err_u_hist(4)
        err_u_hist(4) = err_u
        if (t_step > 3) then
            d_err_u = (err_u_hist(4) - err_u_hist(3))/dt
            ! d_err_u = (1.5_wp*err_u_hist(4) - 2._wp*err_u_hist(3) + 0.5_wp*err_u_hist(2))/dt
            ! d_err_u = (11._wp*err_u_hist(4) - 18._wp*err_u_hist(3) + 9._wp*err_u_hist(2) - 2._wp*err_u_hist(1))/(6._wp*dt)
        else
            d_err_u = 0._wp
        end if

        if (t_step > particle_control_start) then
            particle_bf = particle_bf + cntrl_p%K_Pg*err_u + cntrl_p%K_Dg*d_err_u
            P_inf_ref = P_inf_ref + cntrl_p%K_Pp*err_M
        end if

        $:GPU_UPDATE(device='[particle_bf, P_inf_ref]')

        if (forcing_wrt .and. proc_rank == 0) then
            print *, 'CONTROL:', particle_bf, P_inf_ref, rho_avg*u_rel*Dp/mu, Mach, rho_avg, u_avg, cs_avg, Vp_avg
            write (103) particle_bf, P_inf_ref, rho_avg*u_rel*Dp/mu, Mach
            flush (103)
        end if

    end subroutine s_update_controllers

    subroutine s_finalize_additional_forcing_module

        integer :: i

        do i = 1, num_dims + 2
            @:DEALLOCATE(q_periodic_force(i)%sf)
        end do
        @:DEALLOCATE(q_periodic_force)

        @:DEALLOCATE(rho_window)
        @:DEALLOCATE(rhou_window)
        @:DEALLOCATE(rhoe_window)

        if (forcing_wrt .and. proc_rank == 0) then
            close (102)
        end if

        if (forcing_wrt .and. proc_rank == 0) then
            close (103)
        end if

    end subroutine s_finalize_additional_forcing_module

end module m_additional_forcing
