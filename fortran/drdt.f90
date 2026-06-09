subroutine calc_drdt(nlayer, iaer, icld, co2, ch4, n2o, ps, ts, &
        albedo_lw, plev, t, q, o3, cldfrac, cldlwc,cldiwc, tauaer_lw, lw_base, &
        fdl_base, ful_base, drdt)

    ! Planck Jacobian
    ! Dimension: (nlayer+1)x(nlayer+1)
    ! Row: R change; Column: T change
    ! Atm: convergence of LW fluxes; Surf: LW towards the surface
    
    use parkind, only : im=>kind_im, rb=>kind_rb
    use parrrtm, only:nbndlw
    use parrrsw, only:jpband, jpb1, jpb2
    
    implicit none

    integer(kind=im), intent(in) :: nlayer   ! number of layers
    integer(kind=im), intent(in) :: iaer     ! 0-no aerosol, 10-aerosol
    integer(kind=im), intent(in) :: icld     ! 0-no cloud, 2-cloud
    real(kind=rb), intent(in)    :: co2          ! co2 ppmv
    real(kind=rb), intent(in)    :: ch4          ! ch4 ppmv
    real(kind=rb), intent(in)    :: n2o          ! n2o ppmv
    real(kind=rb), intent(in)    :: ps           ! surface pressure (hPa)
    real(kind=rb), intent(in)    :: ts           ! surface temperature
    real(kind=rb), intent(in)    :: albedo_lw    ! surface albedo for lw
    real(kind=rb), intent(in)    :: plev(nlayer) ! pressures of ERA
    real(kind=rb), intent(in)    :: t(nlayer)    ! atmospheric temperature
    real(kind=rb), intent(in)    :: q(nlayer)    ! specific humidity (kg/kg)
    real(kind=rb), intent(in)    :: o3(nlayer)   ! ozone ppmv
    real(kind=rb), intent(in)    :: cldfrac(nlayer)! cloud fraction (0-1)
    real(kind=rb), intent(in)    :: cldlwc(nlayer) ! cloud liquid water path (kg/kg)
    real(kind=rb), intent(in)    :: cldiwc(nlayer) ! cloud ice water path (kg/kg)
    real(kind=rb), intent(in)    :: tauaer_lw(nlayer,nbndlw)    ! optical depth for longwave each band
    real(kind=rb), intent(in)    :: lw_base(nlayer)
    real(kind=rb), intent(in)    :: fdl_base(nlayer+1)
    real(kind=rb), intent(in)    :: ful_base(nlayer+1)
    real(kind=rb), intent(out)   :: drdt(nlayer+1,nlayer+1)

    integer(kind=im):: ilayer
    real(kind=rb)   :: fdl_p(nlayer+1), fdl_m(nlayer+1)
    real(kind=rb)   :: ful_p(nlayer+1), ful_m(nlayer+1)
    real(kind=rb)   :: htr_lw_1k(nlayer)
    real(kind=rb)   :: lw_p(nlayer), lw_m(nlayer)
    real(kind=rb)   :: t_pert(nlayer), ts_p, ts_m
    real(kind=rb), parameter :: dT_half = 0.5
    logical :: use_centered_drdt

    ! Centered finite-difference Planck probe (T_j ± 0.5K instead of T_j + 1K).
    ! Cancels the leading R_TT term in the FD expansion:
    !   J^onesided = R_T|_base + 0.5·R_TT|_base + O(1)
    !   J^centered = R_T|_base + (1/24)·R_TTT|_base + O(...)
    ! Cost: 2× rad_driver_lw calls in this subroutine (vs 1× for one-sided).
    inquire(file='data_prep/drdt_centered.flag', exist=use_centered_drdt)

    if (use_centered_drdt) then
        ! Atmospheric T_j ± 0.5K
        do ilayer = 1, nlayer
            ! +0.5K
            t_pert(1:nlayer) = t(1:nlayer)
            t_pert(ilayer) = t(ilayer) + dT_half
            call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts, albedo_lw, &
                plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
                fdl_p, ful_p, htr_lw_1k, lw_p)
            ! -0.5K
            t_pert(ilayer) = t(ilayer) - dT_half
            call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts, albedo_lw, &
                plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
                fdl_m, ful_m, htr_lw_1k, lw_m)
            ! Centered FD over total span 1.0K
            drdt(ilayer, 1:nlayer) = lw_p(1:nlayer) - lw_m(1:nlayer)
            drdt(ilayer, nlayer+1) = (fdl_p(nlayer+1)-ful_p(nlayer+1)) &
                                   - (fdl_m(nlayer+1)-ful_m(nlayer+1))
        end do
        ! Surface ts ± 0.5K — perturb ONLY ts (NOT t_atm[nlayer]).
        ! RRTMG's rad_driver_lw uses tint_cfram(nlayer)=ts directly at the
        ! surface interface (line ~90 of rad_driver_lw.f90), so ts and
        ! t_atm[nlayer] are independent inputs. A previous "Fix #1" bundled
        ! them, which contaminated drdt's surface row → wrong sign in
        ! dT_lhflx/dT_shflx at lev=1000hPa (atm bottom) when surface forcing
        ! propagates through J^{-1}. Reverted 2026-05-12 to preserve the
        ! physical Planck Jacobian: ∂R/∂T_s and ∂R/∂T_atm[nlayer] stay separate.
        t_pert(1:nlayer) = t(1:nlayer)
        ts_p = ts + dT_half
        call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts_p, albedo_lw, &
            plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
            fdl_p, ful_p, htr_lw_1k, lw_p)
        ts_m = ts - dT_half
        call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts_m, albedo_lw, &
            plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
            fdl_m, ful_m, htr_lw_1k, lw_m)
        drdt(nlayer+1, 1:nlayer) = lw_p(1:nlayer) - lw_m(1:nlayer)
        drdt(nlayer+1, nlayer+1) = (fdl_p(nlayer+1)-ful_p(nlayer+1)) &
                                 - (fdl_m(nlayer+1)-ful_m(nlayer+1))
    else
        ! Default: one-sided FD (T + 1K), subtract base
        do ilayer = 1, nlayer
            t_pert(1:nlayer) = t(1:nlayer)
            t_pert(ilayer) = t(ilayer) + 1.0
            call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts, albedo_lw, &
                plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
                fdl_p, ful_p, htr_lw_1k, lw_p)
            drdt(ilayer, 1:nlayer) = lw_p(1:nlayer) - lw_base(1:nlayer)
            drdt(ilayer, nlayer+1) = fdl_p(nlayer+1) - ful_p(nlayer+1) - (fdl_base(nlayer+1) - ful_base(nlayer+1))
        end do
        ! Surface ts +1K — perturb ONLY ts (NOT t_atm[nlayer]). See note
        ! above in centered_drdt branch for rationale (reverted Fix #1).
        t_pert(1:nlayer) = t(1:nlayer)
        ts_p = ts + 1.0
        call rad_driver_lw(nlayer, iaer, icld, co2, ch4, n2o, ps, ts_p, albedo_lw, &
            plev, t_pert, q, o3, cldfrac, cldlwc, cldiwc, tauaer_lw, &
            fdl_p, ful_p, htr_lw_1k, lw_p)
        drdt(nlayer+1, 1:nlayer) = lw_p(1:nlayer) - lw_base(1:nlayer)
        drdt(nlayer+1, nlayer+1) = fdl_p(nlayer+1) - ful_p(nlayer+1) - (fdl_base(nlayer+1) - ful_base(nlayer+1))
    end if

end subroutine
     


