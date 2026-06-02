program cfram_fu
!  pyCFRAM single-column driver using Fu radiation (apple-to-apple with
!  OLD CFRAM raw/CFRAM.zip). Mirrors cfram_rrtmg.f90 I/O contract:
!    - Reads same data_prep/*.dat input files (plev/t/q/o3/cc/clwc/ciwc/...)
!    - Writes same data_output/frc_*.dat + drdt_inv.dat (same record format)
!    - Skips aerosol perturbation entirely (no Δaerosol partial)
!
!  Fu radiation interface (cas_fu_radiation.f + fu_helpers.f):
!    - State passed via F77 common blocks: /atmosp/, /clouds/, /umcon/, /zdim/
!    - S_R_cloudy(u0, as, ss, pts, rad, area, sw, lw, water, ice, iseed, no_cloud_out)
!      = full-sky radiation (Monte Carlo cloud sub-sampling, 100 sub-cols/call)
!    - Cloud effective radii hardcoded: pre=15 μm (water), pde=40 μm (ice)
!
!  Runtime nlev: derived from data_prep/plev.dat file size (single binary).

  implicit none
  ! Inline para.file constants (avoid F77 include in free-form)
  integer, parameter :: mbs=6, mbir=12, mb=18, nc=8
  integer, parameter :: nlat=1, nlon=1, ntime=1
  real, parameter    :: scon_const=1360.98

  integer :: nlev, plev_fsize
  integer :: k, kk, ipert, iseed, ib
  integer :: no_cloud_out(100,100)
  ! Per-case nv tracking (apple-to-apple OLD CFRAM GW-cfram.f L268-289):
  ! all partial cases (0,2,3,4,5,6,7,8) use nv_base_used; case 9 (full warm)
  ! uses nv_warm_used. Required for correct frc index alignment when
  ! nv_warm ≠ nv_base.
  integer :: nv_base_used, nv_warm_used, nv_overlap
  integer :: nv1_base, nv1_warm

  real, allocatable :: plev(:)
  real :: ts_base, ts_warm, ts_use, ts_use_p
  real :: ps_base, ps_warm, ps_hpa
  real :: huss_base, huss_warm   ! 2m specific humidity (CMIP6 huss); apple-
                                 ! to-apple OLD CFRAM raw/CFRAM.zip ph(nv1).
                                 ! Sentinel |huss| > 900 → fall back to ph(nv).
  real :: solin_base, solin_warm
  real :: ssrd_base_v, ssrd_warm_v, ssru_base_v, ssru_warm_v
  real :: albedo_base_v, albedo_warm_v
  real :: co2ppmv_base, co2ppmv_warm, ch4ppmv, n2oppmv

  real, allocatable :: t_base_v(:), t_warm_v(:), q_base_v(:), q_warm_v(:)
  real, allocatable :: o3_base_v(:), o3_warm_v(:)
  real, allocatable :: cc_base_v(:), cc_warm_v(:)
  real, allocatable :: clwc_base_v(:), clwc_warm_v(:)
  real, allocatable :: ciwc_base_v(:), ciwc_warm_v(:)

  ! Per-state radiation
  real, allocatable :: rad_base(:), lw_base(:), sw_base(:)
  real, allocatable :: rad_warm(:), lw_warm(:), sw_warm(:)
  real, allocatable :: rad_pert(:), lw_pert(:), sw_pert(:)
  real, allocatable :: rad_full(:), lw_full(:), sw_full(:)
  real, allocatable :: rad_co2(:),  rad_q(:),    rad_o3(:)
  real, allocatable :: rad_solar(:),rad_albedo(:),rad_cloud(:)
  real, allocatable :: lw_cloud(:), sw_cloud(:)
  real, allocatable :: rad_ts(:)     ! case 5: ts-only perturbation (atm at base, sfc T at warm)

  ! Forcing outputs (length nv1 = nlev+1; index nv1 = surface)
  real, allocatable :: frc_warm_o(:), frc_co2_o(:), frc_q_o(:)
  real, allocatable :: frc_o3_o(:),  frc_solar_o(:), frc_albedo_o(:)
  real, allocatable :: frc_cloud_o(:), frc_full_o(:), frc_ts_o(:)
  real, allocatable :: frc_cloud_lw_o(:), frc_cloud_sw_o(:)
  real, allocatable :: frc_zero(:)   ! for frc_aerosol (unused but required)

  ! Planck matrix
  real, allocatable :: drdt(:,:), drdt_inv(:,:)
  real, allocatable :: rad_t_pert(:)
  ! LAPACK args MUST be INTEGER*4 (LP64 MKL/LAPACK), not the default INTEGER*8
  ! produced by -fdefault-integer-8 / -i8. Mismatch crashes DGETRF.
  integer(kind=4), allocatable :: ipiv(:)
  integer(kind=4) :: info, n_lapack
  real, allocatable :: work(:)

  ! Fu common blocks (must align with cas_fu_radiation.f / fu_helpers.f)
  integer :: nv1, nv, ndfs, mdfs, ndfs4, ndfs2
  common /zdim/ nv1, nv, ndfs, mdfs, ndfs4, ndfs2
  real :: pp(100), pt(100), ptc(100), ph(100), po(100)
  common /atmosp/ pp, pt, ptc, ph, po
  real :: pre(100), plwc(100), pde(100), piwc(100)
  common /clouds/ pre, plwc, pde, piwc
  real :: prwc(100)
  common /rains/ prwc
  real :: pgwc(100)
  common /graups/ pgwc
  real :: umco2, umch4, umn2o
  common /umcon/ umco2, umch4, umn2o
  real :: fds(100), fus(100), dts(100), fdir(100), fuir(100), dtir(100)
  real :: fd(100), fu(100), ht(100)
  common /radiat/ fds, fus, dts, fdir, fuir, dtir, fd, fu, ht

  ! Apple-to-apple OLD CFRAM: TWO no_cloud MC patterns are needed because the
  ! sub-column overlap realization must be consistent with the cloud fraction
  ! profile used by each case. OLD CFRAM:
  !   GW-base/GW-co2/GW-wv/GW-o3/GW-solar/GW-albedo/GW-ta/GW-ts/GW-drdt
  !     → open base_no_cloud_out_*.dat (base cc consistent)
  !   GW-warm/GW-cloud
  !     → open warm_no_cloud_out_*.dat (warm cc consistent)
  ! We mirror this in-memory: case 0 generates base from cc_base, case 1
  ! generates warm from cc_warm; cases 8/9 use warm; everything else uses base.
  ! The /cloud_fix/ common block is shared with fu_helpers.f (S_R_cloudy reads
  ! no_cloud_fixed when use_fixed_cloud==1). The two pattern buffers below
  ! are LOCAL to the main program and copied into no_cloud_fixed per case.
  integer :: use_fixed_cloud, no_cloud_fixed(100,100)
  common /cloud_fix/ use_fixed_cloud, no_cloud_fixed
  integer :: no_cloud_base_pat(100,100), no_cloud_warm_pat(100,100)
  integer :: have_base_pat, have_warm_pat

  real :: as(mbs)
  real :: area_col(100), water_col(100), ice_col(100)
  real :: pmean, dry_density
  real :: sol_use, alb_use

  ch4ppmv = 1.6
  n2oppmv = 0.28
  iseed = -12345   ! ran3 re-inits internal state when idum < 0
  use_fixed_cloud = 0   ! case 0 generates; flipped to 1 after
  have_base_pat = 0     ! base MC pattern (cc_base-consistent) not loaded yet
  have_warm_pat = 0     ! warm MC pattern (cc_warm-consistent) not loaded yet

  !---------------- Determine nlev from plev.dat size ----------
  inquire(file='data_prep/plev.dat', size=plev_fsize)
  if (plev_fsize <= 0) then
     print *, 'ERROR: data_prep/plev.dat missing'
     stop 1
  end if
  nlev = plev_fsize / 8
  print *, 'cfram_fu_1col: nlev =', nlev

  !---------------- Allocate ----------
  allocate(plev(nlev))
  allocate(t_base_v(nlev), t_warm_v(nlev))
  allocate(q_base_v(nlev), q_warm_v(nlev))
  allocate(o3_base_v(nlev), o3_warm_v(nlev))
  allocate(cc_base_v(nlev), cc_warm_v(nlev))
  allocate(clwc_base_v(nlev), clwc_warm_v(nlev))
  allocate(ciwc_base_v(nlev), ciwc_warm_v(nlev))
  allocate(rad_base(nlev+1),  lw_base(nlev+1),  sw_base(nlev+1))
  allocate(rad_warm(nlev+1),  lw_warm(nlev+1),  sw_warm(nlev+1))
  allocate(rad_pert(nlev+1),  lw_pert(nlev+1),  sw_pert(nlev+1))
  allocate(rad_full(nlev+1),  lw_full(nlev+1),  sw_full(nlev+1))
  allocate(rad_co2(nlev+1),   rad_q(nlev+1),    rad_o3(nlev+1))
  allocate(rad_solar(nlev+1), rad_albedo(nlev+1),rad_cloud(nlev+1))
  allocate(lw_cloud(nlev+1),  sw_cloud(nlev+1))
  allocate(rad_ts(nlev+1))
  allocate(frc_warm_o(nlev+1), frc_co2_o(nlev+1), frc_q_o(nlev+1))
  allocate(frc_o3_o(nlev+1),   frc_solar_o(nlev+1), frc_albedo_o(nlev+1))
  allocate(frc_cloud_o(nlev+1),frc_full_o(nlev+1), frc_ts_o(nlev+1))
  allocate(frc_cloud_lw_o(nlev+1), frc_cloud_sw_o(nlev+1))
  allocate(frc_zero(nlev+1))
  allocate(drdt(nlev+1,nlev+1), drdt_inv(nlev+1,nlev+1))
  allocate(rad_t_pert(nlev+1), ipiv(nlev+1), work(nlev+1))
  frc_zero = 0.0

  !---------------- Read inputs ----------
  open(99, file='data_prep/plev.dat', form='unformatted', access='direct', recl=nlev*8)
  read(99,rec=1) plev
  close(99)

  call read1d('data_prep/t_base.dat',   nlev, t_base_v)
  call read1d('data_prep/t_warm.dat',   nlev, t_warm_v)
  call read1d('data_prep/hus_base.dat', nlev, q_base_v)
  call read1d('data_prep/hus_warm.dat', nlev, q_warm_v)
  call read1d('data_prep/O3_base.dat',  nlev, o3_base_v)
  call read1d('data_prep/O3_warm.dat',  nlev, o3_warm_v)
  call read1d('data_prep/cc_base.dat',  nlev, cc_base_v)
  call read1d('data_prep/cc_warm.dat',  nlev, cc_warm_v)
  call read1d('data_prep/clwc_base.dat',nlev, clwc_base_v)
  call read1d('data_prep/clwc_warm.dat',nlev, clwc_warm_v)
  call read1d('data_prep/ciwc_base.dat',nlev, ciwc_base_v)
  call read1d('data_prep/ciwc_warm.dat',nlev, ciwc_warm_v)
  call read1scalar('data_prep/skt_base.dat', ts_base)
  call read1scalar('data_prep/skt_warm.dat', ts_warm)
  call read1scalar('data_prep/sp_base.dat',  ps_base)
  call read1scalar('data_prep/sp_warm.dat',  ps_warm)
  call read1scalar_optional('data_prep/huss_base.dat', huss_base, -999.0d0)
  call read1scalar_optional('data_prep/huss_warm.dat', huss_warm, -999.0d0)
  call read1scalar('data_prep/solarin_base.dat', solin_base)
  call read1scalar('data_prep/solarin_warm.dat', solin_warm)
  call read1scalar('data_prep/ssrd_base.dat', ssrd_base_v)
  call read1scalar('data_prep/ssrd_warm.dat', ssrd_warm_v)
  call read1scalar('data_prep/ssru_base.dat', ssru_base_v)
  call read1scalar('data_prep/ssru_warm.dat', ssru_warm_v)
  call read1scalar('data_prep/co2_b.dat', co2ppmv_base)
  call read1scalar('data_prep/co2_w.dat', co2ppmv_warm)

  print *, '  ts_base=', ts_base, ' ts_warm=', ts_warm
  print *, '  ps_base=', ps_base, ' Pa  solin_base=', solin_base
  print *, '  CO2 base=', co2ppmv_base, ' warm=', co2ppmv_warm

  ! albedo
  if (ssrd_base_v > 1.0) then
     albedo_base_v = ssru_base_v / ssrd_base_v
  else
     albedo_base_v = 0.3
  end if
  if (ssrd_warm_v > 1.0) then
     albedo_warm_v = ssru_warm_v / ssrd_warm_v
  else
     albedo_warm_v = 0.3
  end if

  !---------------- Setup Fu vertical grid ----------
  ! pp/nv1/nv now set per-call inside set_state (varies with state's ps —
  ! base uses ps_base, warm-state cases (1, 9) use ps_warm).
  ! We just print the base ps for reference here.
  ps_hpa = ps_base / 100.0
  print *, '  ps_base=', ps_hpa, ' hPa  ps_warm=', ps_warm/100.0, ' hPa'

  umco2 = co2ppmv_base
  umch4 = ch4ppmv
  umn2o = n2oppmv

  ! Optional: preload OLD CFRAM's base/warm no_cloud MC patterns (diagnostic
  ! mode for apple-to-apple comparison). For self-runs (cesm2/eh* cases) these
  ! files don't exist → patterns are generated fresh inside the partial-case
  ! loop: case 0 generates base from cc_base, case 1 generates warm from cc_warm.
  call load_no_cloud_seed('data_prep/base_no_cloud_seed.dat', &
       have_base_pat, no_cloud_base_pat)
  call load_no_cloud_seed('data_prep/warm_no_cloud_seed.dat', &
       have_warm_pat, no_cloud_warm_pat)

  !---------------- Per-perturbation radiation calls ----------
  ! Apple-to-apple with raw/CFRAM.zip OLD CFRAM (GW-base.f, GW-co2.f, …):
  ! each case is an independent program in OLD CFRAM. We mirror that fresh-
  ! program semantics by zeroing all Fu common-block arrays before every
  ! set_state call (no cross-case state contamination).
  do ipert = 0, 9
     call reset_fu_state()
     call set_state(ipert, nlev, mbs, &
          t_base_v, t_warm_v, q_base_v, q_warm_v, &
          o3_base_v, o3_warm_v, cc_base_v, cc_warm_v, &
          clwc_base_v, clwc_warm_v, ciwc_base_v, ciwc_warm_v, &
          ts_base, ts_warm, solin_base, solin_warm, &
          albedo_base_v, albedo_warm_v, &
          co2ppmv_base, co2ppmv_warm, plev, ps_base, ps_warm, &
          huss_base, huss_warm, &
          area_col, water_col, ice_col, &
          as, ts_use, sol_use)

     ! Pick MC pattern for this case. Cases with warm cloud state (1=warm,
     ! 8=cloud, 9=full) need warm_no_cloud (cc_warm-consistent); everything
     ! else (0,2,3,4,5,6,7) uses base_no_cloud (cc_base-consistent). If the
     ! pattern hasn't been generated/loaded yet, fall back to ran3
     ! generation inside S_R_cloudy from whatever cc is currently set.
     if (ipert == 1 .or. ipert == 8 .or. ipert == 9) then
        if (have_warm_pat == 1) then
           no_cloud_fixed(:,:) = no_cloud_warm_pat(:,:)
           use_fixed_cloud = 1
        else
           use_fixed_cloud = 0   ! generate fresh from cc_warm
        end if
     else
        if (have_base_pat == 1) then
           no_cloud_fixed(:,:) = no_cloud_base_pat(:,:)
           use_fixed_cloud = 1
        else
           use_fixed_cloud = 0   ! generate fresh from cc_base (case 0 path)
        end if
     end if

     ! CRITICAL: reset iseed to same value before each call so the Monte
     ! Carlo cloud overlap pattern is reproducible. With use_fixed_cloud=1
     ! the pattern is copied from no_cloud_fixed; with =0 it's regenerated
     ! via ran3 from current cc. Then frc_X = R_X − R_base contains only
     ! the X-perturbation signal, not stochastic cloud noise.
     iseed = -12345
     call S_R_cloudy(sol_use/scon_const, as, scon_const, ts_use, &
          rad_pert, area_col, sw_pert, lw_pert, &
          water_col, ice_col, iseed, no_cloud_out)

     ! Capture freshly-generated pattern (only if not already loaded from disk).
     ! Apple-to-apple with OLD CFRAM: GW-base.f writes base_no_cloud_out_*.dat,
     ! GW-warm.f writes warm_no_cloud_out_*.dat for downstream cases to read.
     if (ipert == 0 .and. have_base_pat == 0) then
        no_cloud_base_pat(:,:) = no_cloud_out(:,:)
        have_base_pat = 1
     end if
     if (ipert == 1 .and. have_warm_pat == 0) then
        no_cloud_warm_pat(:,:) = no_cloud_out(:,:)
        have_warm_pat = 1
     end if

     select case (ipert)
     case (0); rad_base = rad_pert; lw_base = lw_pert; sw_base = sw_pert
              nv_base_used = nv;  nv1_base = nv1
     case (1); rad_warm = rad_pert; lw_warm = lw_pert; sw_warm = sw_pert
     case (2); rad_co2  = rad_pert
     case (3); rad_q    = rad_pert
     case (4); rad_o3   = rad_pert
     case (5); rad_ts   = rad_pert    ! sfc T perturbation only (atm at base)
     case (6); rad_solar  = rad_pert
     case (7); rad_albedo = rad_pert
     case (8); rad_cloud = rad_pert; lw_cloud = lw_pert; sw_cloud = sw_pert
     case (9); rad_full  = rad_pert; lw_full = lw_pert; sw_full = sw_pert
              nv_warm_used = nv;  nv1_warm = nv1
     end select
  end do

  !---------------- Compute frc_X = rad_X − rad_base ----------
  ! Layout: frc_X(1..nv)=atm rows, frc_X(nv+1..nlev)=0 padding,
  ! frc_X(nlev+1)=surface (fixed at LAST index for Python compat).
  !
  ! Apple-to-apple OLD CFRAM (GW-cfram.f L268-289): use per-case nv for index
  ! alignment. Partial cases all use nv_base_used (since they used ps_base);
  ! case 9 (full warm) used nv_warm_used (ps_warm), which can differ. The
  ! buggy version used the GLOBAL nv from /zdim/ (= nv_warm after case 9),
  ! producing wrong subtraction at indices nv_overlap+1..nv_warm and a wrong
  ! surface index when nv_warm ≠ nv_base — that's the source of WV/CO2
  ! magnitude under-estimation and DYN noise.
  frc_co2_o = 0.0;     frc_q_o = 0.0;       frc_o3_o = 0.0
  frc_solar_o = 0.0;   frc_albedo_o = 0.0;  frc_cloud_o = 0.0
  frc_warm_o = 0.0;    frc_full_o = 0.0;    frc_ts_o = 0.0
  frc_cloud_lw_o = 0.0; frc_cloud_sw_o = 0.0

  ! --- partial cases: both rad_X and rad_base have nv = nv_base_used ---
  do k = 1, nv_base_used
     frc_co2_o(k)      = rad_co2(k)    - rad_base(k)
     frc_q_o(k)        = rad_q(k)      - rad_base(k)
     frc_o3_o(k)       = rad_o3(k)     - rad_base(k)
     frc_solar_o(k)    = rad_solar(k)  - rad_base(k)
     frc_albedo_o(k)   = rad_albedo(k) - rad_base(k)
     frc_cloud_o(k)    = rad_cloud(k)  - rad_base(k)
     frc_warm_o(k)     = rad_warm(k)   - rad_base(k)
     frc_ts_o(k)       = rad_ts(k)     - rad_base(k)
     frc_cloud_lw_o(k) = lw_cloud(k)   - lw_base(k)
     frc_cloud_sw_o(k) = sw_cloud(k)   - sw_base(k)
  end do
  ! Surface for partial cases: at nv1_base (= nv_base_used+1)
  frc_co2_o(nlev+1)      = rad_co2(nv1_base)    - rad_base(nv1_base)
  frc_q_o(nlev+1)        = rad_q(nv1_base)      - rad_base(nv1_base)
  frc_o3_o(nlev+1)       = rad_o3(nv1_base)     - rad_base(nv1_base)
  frc_solar_o(nlev+1)    = rad_solar(nv1_base)  - rad_base(nv1_base)
  frc_albedo_o(nlev+1)   = rad_albedo(nv1_base) - rad_base(nv1_base)
  frc_cloud_o(nlev+1)    = rad_cloud(nv1_base)  - rad_base(nv1_base)
  frc_warm_o(nlev+1)     = rad_warm(nv1_base)   - rad_base(nv1_base)
  frc_ts_o(nlev+1)       = rad_ts(nv1_base)     - rad_base(nv1_base)
  frc_cloud_lw_o(nlev+1) = lw_cloud(nv1_base)   - lw_base(nv1_base)
  frc_cloud_sw_o(nlev+1) = sw_cloud(nv1_base)   - sw_base(nv1_base)

  ! --- case 9 (full warm): rad_full has nv = nv_warm_used, possibly ≠ nv_base ---
  ! Atmospheric rows: only the overlap region is physically meaningful.
  ! OLD CFRAM GW-cfram.f L268-289 uses min(nv_base, nv_warm) → exit on -999.
  nv_overlap = min(nv_base_used, nv_warm_used)
  do k = 1, nv_overlap
     frc_full_o(k) = rad_full(k) - rad_base(k)
  end do
  ! Surface: each case at its own nv1. Both represent the surface row, even if
  ! ps_warm ≠ ps_base (slightly different surface pressures, but conceptually
  ! both are "surface").
  frc_full_o(nlev+1) = rad_full(nv1_warm) - rad_base(nv1_base)

  !---------------- Compute drdt (T+1K each layer at base) ----------
  print *, '  Computing drdt: ', nv1, ' perturbations of +1K each layer'
  ! Apple-to-apple with OLD CFRAM GW-drdt.f: drdt computed on a fresh base
  ! state (clean common blocks, ps_base). We MUST also re-establish lw_base
  ! by re-calling S_R_cloudy on this fresh base — otherwise the lw_base
  ! saved from the loop above used a state with non-zero common-block tail
  ! left from initial program startup.
  call reset_fu_state()
  call set_state(0, nlev, mbs, &
       t_base_v, t_warm_v, q_base_v, q_warm_v, &
       o3_base_v, o3_warm_v, cc_base_v, cc_warm_v, &
       clwc_base_v, clwc_warm_v, ciwc_base_v, ciwc_warm_v, &
       ts_base, ts_warm, solin_base, solin_warm, &
       albedo_base_v, albedo_warm_v, &
       co2ppmv_base, co2ppmv_warm, plev, ps_base, ps_warm, &
       huss_base, huss_warm, &
       area_col, water_col, ice_col, &
       as, ts_use, sol_use)
  ! drdt is computed at base state — apple-to-apple with OLD GW-drdt.f, which
  ! reads base_no_cloud_out. Force base MC pattern for all 19 S_R_cloudy calls
  ! in this block (1 base + 17 atm-T perturbations + 1 surface-T perturbation).
  no_cloud_fixed(:,:) = no_cloud_base_pat(:,:)
  use_fixed_cloud = 1
  iseed = -12345
  call S_R_cloudy(sol_use/scon_const, as, scon_const, ts_use, &
       rad_base, area_col, sw_base, lw_base, &
       water_col, ice_col, iseed, no_cloud_out)

  ! atm layers: perturb pt(k) by +1, recompute (cloud sky), drdt[:,k] = LW response.
  ! Reset iseed before each S_R_cloudy so identical MC cloud overlap → drdt
  ! contains only T-perturbation signal, not stochastic cloud noise.
  do k = 1, nv
     pt(k) = pt(k) + 1.0
     iseed = -12345
     call S_R_cloudy(sol_use/scon_const, as, scon_const, ts_use, &
          rad_t_pert, area_col, sw_pert, lw_pert, &
          water_col, ice_col, iseed, no_cloud_out)
     pt(k) = pt(k) - 1.0
     do kk = 1, nv1
        drdt(kk, k) = lw_pert(kk) - lw_base(kk)
     end do
  end do
  ! Surface T perturbation — perturb ONLY ts (the skin / surface emission T,
  ! passed as `pts` argument). pt(nv1) is the atm-bottom interface T (= pt(nv)
  ! per set_state below) and is NOT physically tied to ts skin. This decouples
  ! skin T and air bottom T like RRTMG does — OLD CFRAM's bundling
  ! (perturbing pt(nv1) and pts together when k0=nv1) was a side effect of
  ! convention `pt(nv1) := ts`, not an intentional physical model. 2026-05-12
  ! reverted to pure ts perturbation per user direction.
  ts_use_p = ts_use + 1.0
  ! pt(nv1) unchanged — pure ts (skin only) perturbation
  iseed = -12345
  call S_R_cloudy(sol_use/scon_const, as, scon_const, ts_use_p, &
       rad_t_pert, area_col, sw_pert, lw_pert, &
       water_col, ice_col, iseed, no_cloud_out)
  do kk = 1, nv1
     drdt(kk, nv1) = lw_pert(kk) - lw_base(kk)
  end do

  !---------------- Diagnostic: dump drdt before inversion ----------
  ! drdt(kk, k) = lw_pert(kk) - lw_base(kk) when pt(k) is +1K (LW only).
  ! For surface-row diagonal: drdt(nv1, nv1) ≈ -4σTs³ (Stefan-Boltzmann),
  ! since +1K skin T → ULW emission +4σTs³ W/m² → net surface LW flux down by
  ! that amount. With Fix #2 also perturbing pt(nv1), atm bottom DLW partially
  ! offsets ULW increase. Expected magnitude ~ -5 to -3 W/m²/K.
  open(97, file='data_output/drdt.dat', access='sequential', &
       form='unformatted', status='replace')
  write(97) int(nv, kind=8)
  write(97) ((drdt(k, kk), kk=1, nv1), k=1, nv1)
  close(97)

  !---------------- Invert drdt (LAPACK) ----------
  ! drdt_inv is allocated (nlev+1, nlev+1) but only the upper-left (nv1, nv1)
  ! block is meaningful (when ps low, sub-surface levels excluded). LAPACK LDA
  ! must be the ACTUAL leading dimension (nlev+1), not the matrix size (nv1).
  drdt_inv = drdt
  n_lapack = int(nv1, kind=4)
  call dgetrf(n_lapack, n_lapack, drdt_inv, int(nlev+1, kind=4), ipiv, info)
  if (info /= 0) then
     print *, 'ERROR: dgetrf info=', info; stop 2
  end if
  call dgetri(n_lapack, drdt_inv, int(nlev+1, kind=4), ipiv, work, int(nlev+1, kind=4), info)
  if (info /= 0) then
     print *, 'ERROR: dgetri info=', info; stop 3
  end if

  !---------------- Write outputs ----------
  ! frc_X arrays are length nlev+1 (atm 1..nv + zero padding nv+1..nlev + surface at nlev+1).
  ! Python's run_parallel_python.py reads ALL nlev+1 values and treats index NLEV (last)
  ! as surface row. drdt_inv is dynamic size (nv+1)×(nv+1), header indicates nv.
  call ensure_dir('data_output')
  call write_frc('data_output/frc_co2.dat',      frc_co2_o,      nlev+1)
  call write_frc('data_output/frc_q.dat',        frc_q_o,        nlev+1)
  call write_frc('data_output/frc_o3.dat',       frc_o3_o,       nlev+1)
  call write_frc('data_output/frc_solar.dat',    frc_solar_o,    nlev+1)
  call write_frc('data_output/frc_albedo.dat',   frc_albedo_o,   nlev+1)
  call write_frc('data_output/frc_cloud.dat',    frc_cloud_o,    nlev+1)
  call write_frc('data_output/frc_aerosol.dat',  frc_zero,       nlev+1)
  call write_frc('data_output/frc_warm.dat',     frc_warm_o,     nlev+1)
  call write_frc('data_output/frc_full.dat',     frc_full_o,     nlev+1)
  call write_frc('data_output/frc_cloud_lw.dat', frc_cloud_lw_o, nlev+1)
  call write_frc('data_output/frc_cloud_sw.dat', frc_cloud_sw_o, nlev+1)
  call write_frc('data_output/frc_ts.dat',       frc_ts_o,       nlev+1)

  ! drdt_inv: header (nlayer=nv:int8) + body ((nv+1)*(nv+1) doubles, F-order)
  open(98, file='data_output/drdt_inv.dat', access='sequential', &
       form='unformatted', status='replace')
  write(98) int(nv, kind=8)            ! nlayer = nv (effective above-surface count)
  write(98) ((drdt_inv(k, kk), kk=1, nv1), k=1, nv1)
  close(98)

  print *, 'cfram_fu_1col: done'

contains

  subroutine read1d(fname, n, arr)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n
    real, intent(out) :: arr(n)
    open(80, file=fname, form='unformatted', access='direct', recl=n*8)
    read(80, rec=1) arr
    close(80)
  end subroutine

  subroutine read1scalar(fname, x)
    character(len=*), intent(in) :: fname
    real, intent(out) :: x
    real :: buf(1)
    open(80, file=fname, form='unformatted', access='direct', recl=8)
    read(80, rec=1) buf
    close(80)
    x = buf(1)
  end subroutine

  subroutine load_no_cloud_seed(fname, flag, ncf)
    ! Optionally load OLD CFRAM's base_no_cloud_out_*.dat record for one cell
    ! (100×100 INTEGER*4, default ifort = 4-byte int). Cast to our INTEGER*8
    ! (-i8 build). Sets flag=1 if loaded, else 0.
    character(len=*), intent(in)  :: fname
    integer, intent(out) :: flag
    integer, intent(out) :: ncf(100, 100)
    logical :: ex
    integer(kind=4) :: tmp(100, 100)
    integer :: i, j
    inquire(file=fname, exist=ex)
    if (.not. ex) then
       flag = 0
       return
    end if
    open(82, file=fname, form='unformatted', access='direct', recl=100*100*4)
    read(82, rec=1) tmp
    close(82)
    do j = 1, 100
       do i = 1, 100
          ncf(i, j) = int(tmp(i, j), kind=8)
       end do
    end do
    flag = 1
    print *, '  Loaded MC seed: ', trim(fname)
  end subroutine

  subroutine read1scalar_optional(fname, x, default_val)
    ! Like read1scalar but returns default_val if file is missing. Used for
    ! optional inputs (huss) that legacy cases (collaborator's cesm2_4xco2)
    ! don't provide; downstream code checks |x| > 900 sentinel.
    character(len=*), intent(in) :: fname
    real, intent(out) :: x
    real, intent(in)  :: default_val
    logical :: ex
    inquire(file=fname, exist=ex)
    if (ex) then
       call read1scalar(fname, x)
    else
       x = default_val
    end if
  end subroutine

  subroutine write_frc(fname, arr, n)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n
    real, intent(in) :: arr(n)
    open(81, file=fname, form='unformatted', status='replace')
    write(81) arr
    close(81)
  end subroutine

  subroutine ensure_dir(d)
    character(len=*), intent(in) :: d
    logical :: ex
    inquire(file=d, exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p '//trim(d))
  end subroutine

end program cfram_fu


subroutine set_state(ipert, nlev, mbs, &
     tb, tw, qb, qw, o3b, o3w, ccb, ccw, &
     clwb, clww, ciwb, ciww, &
     tsb, tsw, solb, solw, alb_b, alb_w, &
     co2b, co2w, plev, psb, psw, &
     huss_b, huss_w, &
     area_col, water_col, ice_col, &
     as, ts_use, sol_use)
  ! Set Fu common-block state for a given perturbation code
  ! ipert: 0=base, 1=warm-lump, 2=co2, 3=q, 4=o3, 5=ts, 6=solar, 7=albedo, 8=cloud, 9=full
  implicit none
  integer, intent(in) :: ipert, nlev, mbs
  real, intent(in)  :: tb(nlev), tw(nlev), qb(nlev), qw(nlev)
  real, intent(in)  :: o3b(nlev), o3w(nlev), ccb(nlev), ccw(nlev)
  real, intent(in)  :: clwb(nlev), clww(nlev), ciwb(nlev), ciww(nlev)
  real, intent(in)  :: tsb, tsw, solb, solw, alb_b, alb_w
  real, intent(in)  :: co2b, co2w
  real, intent(in)  :: plev(nlev)
  real, intent(in)  :: psb, psw      ! ps_base, ps_warm in Pa
  real, intent(in)  :: huss_b, huss_w  ! 2m specific humidity (sentinel ≤ -900 = absent)
  real, intent(out) :: area_col(100), water_col(100), ice_col(100)
  real, intent(out) :: as(mbs), ts_use, sol_use

  integer :: l, ib, k
  real :: pmean, dry_density, alb_use, ps_use, huss_use

  integer :: nv1, nv, ndfs, mdfs, ndfs4, ndfs2
  common /zdim/ nv1, nv, ndfs, mdfs, ndfs4, ndfs2
  real :: pp(100), pt(100), ptc(100), ph(100), po(100)
  common /atmosp/ pp, pt, ptc, ph, po
  real :: pre(100), plwc(100), pde(100), piwc(100)
  common /clouds/ pre, plwc, pde, piwc
  real :: umco2, umch4, umn2o
  common /umcon/ umco2, umch4, umn2o

  ! Step 1: pick this state's ps. Apple-to-apple with raw/CFRAM.zip OLD CFRAM:
  !   only GW-warm.f reads sp_warm.dat (case 9 = full warm); all other
  !   per-perturbation programs (GW-base, GW-co2, GW-o3, GW-solar, GW-albedo,
  !   GW-cloud, GW-wv, GW-ts) read sp_base.dat.
  ! Case 1 (warm-lump, our addition not in OLD CFRAM) keeps T_base, so ps_base
  ! is the natural choice (T-driven hydrostatic adjustment of ps requires T_warm).
  if (ipert == 9) then
     ps_use = psw / 100.0
  else
     ps_use = psb / 100.0
  end if

  ! Apple-to-apple OLD CFRAM huss mapping:
  !   GW-base / GW-co2 / GW-o3 / GW-ts / GW-solar / GW-albedo / GW-cloud
  !     read huss_base.dat → ph(nv1) = huss_base
  !   GW-wv (case 3, q perturbation):     huss_warm
  !   GW-warm (case 9, full warm):        huss_warm
  ! Case 1 (warm-lump, others_warm + T_base): mirror GW-wv (q is warm).
  if (ipert == 1 .or. ipert == 3 .or. ipert == 9) then
     huss_use = huss_w
  else
     huss_use = huss_b
  end if
  nv = 0
  do k = 1, nlev
     if (plev(k) < ps_use) then
        nv = k
        pp(k) = plev(k)
     else
        exit
     end if
  end do
  nv1 = nv + 1
  pp(nv1) = ps_use
  ndfs = nv;  mdfs = nv1;  ndfs4 = 4*ndfs;  ndfs2 = 2*ndfs

  ! Step 2: fill state arrays for first nv levels per perturbation
  select case (ipert)
  case (0)   ! base
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_b
     umco2 = co2b
  case (1)   ! warm (T base, others warm) — OLD lump
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qw(1:nv);  po(1:nv) = o3w(1:nv)
     area_col(1:nv) = ccw(1:nv);  water_col(1:nv) = clww(1:nv);  ice_col(1:nv) = ciww(1:nv)
     ts_use = tsw;  sol_use = solw;  alb_use = alb_w
     umco2 = co2w
  case (2)   ! co2 only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_b
     umco2 = co2w
  case (3)   ! q only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qw(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_b
     umco2 = co2b
  case (4)   ! o3 only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3w(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_b
     umco2 = co2b
  case (5)   ! ts only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsw;  sol_use = solb;  alb_use = alb_b
     umco2 = co2b
  case (6)   ! solar only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solw;  alb_use = alb_b
     umco2 = co2b
  case (7)   ! albedo only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccb(1:nv);  water_col(1:nv) = clwb(1:nv);  ice_col(1:nv) = ciwb(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_w
     umco2 = co2b
  case (8)   ! cloud only
     pt(1:nv) = tb(1:nv);  ph(1:nv) = qb(1:nv);  po(1:nv) = o3b(1:nv)
     area_col(1:nv) = ccw(1:nv);  water_col(1:nv) = clww(1:nv);  ice_col(1:nv) = ciww(1:nv)
     ts_use = tsb;  sol_use = solb;  alb_use = alb_b
     umco2 = co2b
  case (9)   ! full warm (everything including T AND ps)
     pt(1:nv) = tw(1:nv);  ph(1:nv) = qw(1:nv);  po(1:nv) = o3w(1:nv)
     area_col(1:nv) = ccw(1:nv);  water_col(1:nv) = clww(1:nv);  ice_col(1:nv) = ciww(1:nv)
     ts_use = tsw;  sol_use = solw;  alb_use = alb_w
     umco2 = co2w
  end select

  ! Surface row (l = nv1 = nv+1). ph(nv1) uses CMIP6 huss when available
  ! (apple-to-apple OLD CFRAM raw/CFRAM.zip GW-base.f L322-330); otherwise
  ! falls back to ph(nv) HOLD (legacy cesm2_4xco2 / collaborator data which
  ! has no huss field in surf NC).
  !
  ! pt(nv1) is the AIR temperature at the surface interface — physically the
  ! atm-bottom layer T extrapolated to ps, NOT the skin (ts). OLD CFRAM bundled
  ! pt(nv1) := ts as a convention shortcut, but that conflates skin T with air
  ! T and breaks the Planck Jacobian's surface-row physics (J⁻¹[atm_bottom,sfc]
  ! becomes too tightly coupled). 2026-05-12: decoupled — pt(nv1) HOLDs the atm
  ! bottom T (= pt(nv)), skin T is passed exclusively via the `pts` argument
  ! to S_R_cloudy. Matches RRTMG's tave/tbound separation.
  pt(nv1) = pt(nv)
  if (abs(huss_use) > 900.0) then
     ph(nv1) = ph(nv)              ! sentinel → fallback
  else
     ph(nv1) = huss_use             ! 2m specific humidity
  end if
  po(nv1) = po(nv)
  do l = 1, nv1
     ptc(l) = pt(l)
  end do

  ! Cloud water/ice content in g/m³ (OLD CFRAM convention, GW-warm.f line 367-376)
  do l = 1, nv
     pre(l)  = 15.0
     pde(l)  = 40.0
     pmean = pp(l+1)
     dry_density = pmean * 100.0 / (287.0 * pt(l))
     plwc(l) = 1000.0 * water_col(l) * dry_density
     piwc(l) = 1000.0 * ice_col(l)   * dry_density
  end do

  do ib = 1, mbs
     as(ib) = alb_use
  end do
end subroutine set_state


subroutine reset_fu_state()
  ! Apple-to-apple with raw/CFRAM.zip OLD CFRAM (GW-base.f L287-294 cloud-array
  ! reset per cell), extended to ALL Fu common-block arrays. OLD CFRAM ran
  ! each perturbation as a separate program → fresh BSS-zeroed common blocks.
  ! Our single-binary multi-case driver must explicitly zero between cases,
  ! otherwise stale values at indices > nv (left by a previous case with
  ! larger nv from a different ps) leak into the current case's RT.
  implicit none
  real :: pp(100), pt(100), ptc(100), ph(100), po(100)
  common /atmosp/ pp, pt, ptc, ph, po
  real :: pre(100), plwc(100), pde(100), piwc(100)
  common /clouds/ pre, plwc, pde, piwc
  real :: prwc(100)
  common /rains/ prwc
  real :: pgwc(100)
  common /graups/ pgwc
  real :: fds(100), fus(100), dts(100)
  real :: fdir(100), fuir(100), dtir(100)
  real :: fd(100), fu(100), ht(100)
  common /radiat/ fds, fus, dts, fdir, fuir, dtir, fd, fu, ht
  real :: fdsljh(100), fusljh(100), fdirljh(100), fuirljh(100)
  real :: fuljh(100), fdljh(100), fsljh(100), firljh(100)
  common /radnew/ fdsljh, fusljh, fdirljh, fuirljh, &
                  fuljh, fdljh, fsljh, firljh
  real :: fuir_clr(100), fdir_clr(100), fus_clr(100), fds_clr(100)
  common /clearsky/ fuir_clr, fdir_clr, fus_clr, fds_clr
  real :: fuir_tot(100), fdir_tot(100), fus_tot(100), fds_tot(100)
  common /totalsky/ fuir_tot, fdir_tot, fus_tot, fds_tot

  pp = 0.0; pt = 0.0; ptc = 0.0; ph = 0.0; po = 0.0
  pre = 0.0; plwc = 0.0; pde = 0.0; piwc = 0.0
  prwc = 0.0; pgwc = 0.0
  fds = 0.0; fus = 0.0; dts = 0.0
  fdir = 0.0; fuir = 0.0; dtir = 0.0
  fd = 0.0; fu = 0.0; ht = 0.0
  fdsljh = 0.0; fusljh = 0.0; fdirljh = 0.0; fuirljh = 0.0
  fuljh = 0.0; fdljh = 0.0; fsljh = 0.0; firljh = 0.0
  fuir_clr = 0.0; fdir_clr = 0.0; fus_clr = 0.0; fds_clr = 0.0
  fuir_tot = 0.0; fdir_tot = 0.0; fus_tot = 0.0; fds_tot = 0.0
end subroutine reset_fu_state
