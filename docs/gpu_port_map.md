# MOM6 GPU Port Map

**Scope:** modules currently touched by `!$omp target` / `do concurrent`. Focused on the split-RK2 dynamics core and its direct neighbours. Written to answer: *"if I want to port subroutine X, what's already on the GPU, what isn't, and where are the sync boundaries?"*

**Source of truth:** line numbers cite actual source files at the time of writing. Re-verify before relying on any single directive — the code is moving quickly.

---

## How to read this document

Every array/field gets one of these residency tags:

| Tag | Meaning |
|---|---|
| `[GPU persistent]` | Allocated on device at module init, released at module end. Safe to assume device-resident in any step-level code. |
| `[GPU transient]` | Allocated on device per step (or per sub-step), released at end of scope. Only valid inside the owning region. |
| `[GPU sync]` | Lives on both host and device; kept in sync via `!$omp target update`. Expensive — every access on the "wrong" side forces a transfer. |
| `[CPU only]` | Never mapped. Reading it inside a `target` region is a bug. Typically pointers to other modules' control structures, diagnostic IDs, or logical flags. |
| `[pointer — verify]` | Declared `pointer` or `allocatable` with conditional association — needs a runtime check. |

**Sync-point flags** used in per-module sections:
- **Hot path** — inside the per-timestep loop
- ️ **Wasted roundtrip** — `update from` → trivial CPU work → `update to` on same array. Low-hanging fruit for porting.
- **Gap** — array used in a GPU region without explicit `enter data map`. Either implicitly mapped (fragile) or a latent bug.
- **Confirmed bug** — verified defect, file+line cited.

---

## Big picture: where is data mapped?

The GPU data lifecycle for one MOM6 run looks like this:

```
┌─────────────────────────────────────────────────────────────┐
│ initialize_MOM (MOM.F90 ~2488–3594) │
│ enter data map(to: CS%US, G, G%dx*, G%dy*, G%mask2d*, │
│ G%area*, G%Iarea*, G%bathyT, G%Coriolis*, │
│ GV, GV%Rlay, GV%g_prime, │
│ CS%u, CS%v, CS%h, CS%uh, CS%vh, │
│ CS%uhtr, CS%vhtr, CS%pbv, ...) │
│ enter data map(alloc: CS%eta_av_bc, CS%dyn_split_RK2_CSp,│
│ CS%visc) │
├─────────────────────────────────────────────────────────────┤
│ initialize_dyn_split_RK2 (MOM_dynamics_split_RK2 ~1338–1792)│
│ enter data map(alloc: CS%CAu/v, CS%CAu/v_pred, │
│ CS%PFu/v, CS%diffu/v, │
│ CS%visc_rem_u/v, CS%u/v_accel_bt, │
│ CS%eta, CS%eta_PF, CS%pbce, │
│ CS%uhbt, CS%vhbt, │
│ CS%u_av, CS%v_av, CS%h_av, │
│ CS%continuity_CSp, │
│ CS%PressureForce_CSp, │
│ CS%hor_visc, CS%vertvisc_CSp, │
│ CS%barotropic_CSp) │
│ enter data map(to: CS%taux_bot, CS%tauy_bot) │
├─────────────────────────────────────────────────────────────┤
│ step_MOM_dyn_split_RK2 (~428–1183, called every step) │
│ enter data map(alloc: up, vp, hp, dz, h_tmp, │
│ u_bc_accel, v_bc_accel, eta_pred, │
│ uh_in, vh_in) │
│ … predictor, btstep, corrector … │
│ exit data map(delete/release: same) │
├─────────────────────────────────────────────────────────────┤
│ end_dyn_split_RK2 / MOM_end │
│ exit data map(delete: … most, not all — see §BUGS) │
└─────────────────────────────────────────────────────────────┘
```

**Key property:** `G` (grid), `GV` (vertical grid), and the core state (`u, v, h, uh, vh`) are **persistent** on the device. Everything in the dynamics step assumes they are there. If you port a new routine, you can treat these as free.

**Second key property:** Control structures are mapped *as structs*, not one field at a time, in most cases. But in several places the code maps both the parent struct and its members separately — this is sometimes necessary (allocatables inside a struct don't map transitively), sometimes redundant. The barotropic CS init (lines 6235–6334 of `MOM_barotropic.F90`) does 23 separate per-field maps; consolidation may be worth it.

---

## Module maps

### `src/core/MOM_dynamics_split_RK2.F90`

**Type: `MOM_dyn_split_RK2_CS`** (lines 87–275)

| Field | Shape | Residency | Notes |
|---|---|---|---|
| `CAu, CAv` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1341, exit 2027 |
| `CAu_pred, CAv_pred` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1344, exit 2029 |
| `PFu, PFv` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1347, exit 2031 |
| `PFu_Stokes, PFv_Stokes` | (IBq,J,K) / (I,JBq,K) | [CPU only] | Allocated 1615–1616, never mapped. Only read/written by CPU `Stokes_PGF` at 537, 906. |
| `diffu, diffv` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1338, exit 2025 |
| `visc_rem_u, visc_rem_v` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1608, ** no exit** (deallocated 2037 without `exit data`) |
| `u_accel_bt, v_accel_bt` | (IBq,J,K) / (I,JBq,K) | [GPU persistent] | enter 1614, ** no exit** (dealloc 2036) |
| `eta, eta_PF` | (I,J) | [GPU persistent] | enter 1352/1611, exit 2040 |
| `uhbt, vhbt` | (IBq,J) / (I,JBq) | [GPU persistent] | enter 1605, ** no exit** (dealloc 2035) |
| `pbce` | (I,J,K) | [GPU persistent] | enter 1611, exit 2040 |
| `u_av, v_av, h_av` | (IBq,J,K) / (I,JBq,K) / (I,J,K) | [GPU persistent] | enter 1352, exit 2042 |
| `taux_bot, tauy_bot` | (I,J) pointers | [GPU persistent, read-only copy] | `map(to:)` at 1601, ** never updated and never exited** |
| `hor_visc` | embedded CS | [GPU persistent] | enter 1673, exit 2018 |
| `continuity_CSp` | embedded CS | [GPU persistent] | enter 1657, ** no exit** |
| `CoriolisAdv` | embedded CS | [CPU only] | `CorAdCalc` is called CPU-side (580, 964, 1190) after `update to` of inputs |
| `PressureForce_CSp` | embedded CS | [GPU persistent] | enter 1669, ** no exit** |
| `vertvisc_CSp` | pointer CS | [GPU persistent] | enter 1677, ** no exit** |
| `barotropic_CSp` | embedded CS | [GPU persistent] | enter 1709, exit 2011 |
| `KPP_CSp`, `energetic_PBL_CSp`, `SAL_CSp`, `tides_CSp`, `HA_CSp`, `ALE_CSp` | pointers | [CPU only] | All consumers are CPU-side |
| `OBC` | pointer | [CPU only] | All OBC routines are CPU; triggers frequent `update from/to` pairs |
| `BT_cont`, `diag`, `ADp`, `AD_pred`, `CDp` | pointers | [CPU only] | Diagnostic/pointer scaffolding |

**Transient per-step arrays** (allocated line 428–429, released 1182–1183): `up, vp, hp, dz, h_tmp, u_bc_accel, v_bc_accel, uh_in, vh_in, eta_pred`. Treat as device-only inside `step_MOM_dyn_split_RK2`.

**`do concurrent` kernels:** 21 kernels in `step_MOM_dyn_split_RK2`, lines 431–1179. All live inside the 428/429 target data region. They read/write CS acceleration fields and the transient arrays above, plus intent(inout) args `h, u_inst, v_inst, uh, vh` and grid masks `G%mask2dCu, G%mask2dCv`.

 **Implicit-mapping risk:** `h, u_inst, v_inst, uh, vh, G%mask2d{Cu,Cv}` are touched inside `do concurrent` but never appear in an explicit `enter data`. They work today because they're already mapped upstream in `MOM.F90` (~3008, 3079). If a caller ever invokes `step_MOM_dyn_split_RK2` without the MOM.F90 init having run, this segfaults silently on device. Worth an assertion or explicit `is_device_ptr` clause.

**Sync points (hot path, per timestep):**

| Line | Dir | Vars | Flag | Cause |
|---|---|---|---|---|
| 522 | from | `CS%eta_PF` | | CPU pressure-to-eta conversion |
| 535 | from | `h, u_inst, v_inst, CS%PFu, CS%PFv` | ️ | `Stokes_PGF` is CPU-only |
| 565 | from | `CS%PFu, CS%PFv` | ️ | `open_boundary_zero_normal_flow` is CPU-only |
| 570, 650 | from/to | `eta` | | halo pass (MPI — genuinely unavoidable until halo-on-GPU) |
| 579 | to | `u_av, v_av, h_av, uh, vh` | | `CorAdCalc` is CPU-only (5 large 3D arrays) |
| 597/599 | from/to | `u_bc_accel, v_bc_accel` | ️ | OBC zeroing roundtrip, predictor |
| 797/800 | from/to | `up, vp, h` / `up, vp` | ️ | `vertFPmix` + `vertvisc` CPU-only |
| 848/858 | from/to | `u_av, v_av` | ️ | OBC roundtrip |
| 904 | from | `CS%PFu, CS%PFv, h` | ️ | `Stokes_PGF` (corrector) |
| 984/986 | from/to | `u_bc_accel, v_bc_accel` | ️ | OBC roundtrip (corrector) |
| 1070 | from | `u_inst, v_inst` | ️ | Zero-init of CPU arrays — no data needed, just sync for symmetry |
| 1093/1096 | from/to | `u_inst, v_inst, h` / `u_inst, v_inst` | ️ | `vertFPmix` + `vertvisc` (corrector) |
| 1160 | from | `u_inst, v_inst` | | `radiation_open_bdry_conds` CPU-only |

Non-hot-path transfers also happen at 604–607, 635, 739, 745–750, 940–942, 991–994, 1021, 1033, 1050, 1056–1059 — all guarded by `CS%debug` or a diagnostic ID and only fire if enabled.

---

### `src/core/MOM_barotropic.F90`

**Type: `barotropic_CS`** (lines 116–450+) — large; only GPU-relevant fields shown.

| Field | Shape | Residency | Notes |
|---|---|---|---|
| `frhatu, frhatv` | 3D | [GPU persistent] | enter 867–874 |
| `IDatu, IDatv` | 2D | [GPU persistent] | enter 867–874 |
| `lin_drag_u, lin_drag_v` | 2D | [GPU persistent] | enter 867–874 |
| `ubtav, vbtav` | 2D | [GPU persistent] | enter 867–874 |
| `ubt_IC, vbt_IC` | 2D | [verify] | Allocated in CS but not visible in any explicit `enter data`. Barotropic init region likely covers them; confirm. |
| `eta_cor, eta_cor_bound` | 2D | [GPU persistent] | enter 867–874 |
| `ua_polarity, va_polarity` | 2D | [GPU persistent] | enter 6235–6334 |
| `bathyT, IareaT` | 2D | [GPU persistent] | enter 6235–6334 |
| `dy_Cu, dx_Cv, IdxCu, IdyCv` | 2D | [GPU persistent] | enter 6235–6334 |
| `OBCmask_u, OBCmask_v` | 2D | [GPU persistent] | enter 6235–6334 |
| `D_u_Cor, D_v_Cor, q_D` | 2D | [GPU persistent] | enter 6235–6334 |
| `BT_OBC` | embedded type | [GPU persistent, partial] | Arrays inside `BT_OBC` (Cg_u/v, dZ_u/v, uhbt, vhbt, ubt_outer/vbt_outer, SSH_outer_u/v, u_OBC_type, v_OBC_type) mapped at 867–874. `pass_uv, scalar_pass` (group_pass_type) are [CPU only] — correct. |

**Init region mapping style (lines 6235–6334):** 23 separate `enter data map(to: CS%field)` directives. Functionally correct but verbose. Could consolidate to a handful of structured `map(to: CS%a, CS%b, CS%c, ...)` clauses for readability. Not a bug.

**`do concurrent` kernels:** 100+, concentrated in lines 880–2241. All within the big barotropic stepping target region. Heavy use of `G%mask2dT/Cu/Cv` and `CS%BT_OBC%u_OBC_type/v_OBC_type` conditionals.

**Suspicious sync patterns:**
- Lines 935/1140: `q, DCor_u, DCor_v` from/to roundtrip.
- Lines 1213–1227: nested from/to pairs on `ubt, vbt, uhbt, vhbt, BTCL_u, BTCL_v` inside conditional blocks. ️ candidate.
- Lines 2733–2815: diagnostic dump region, ~40 individual `update from` on PFu/v, Cor_u/v, BT_force_u/v, submerged. Guarded by diagnostic IDs — not hot unless diagnostics enabled.
- ** Lines 1504, 1674–1676:** conditional `update if(...) from(...)` on `uhbt0/vhbt0`, `eta_IC`, `eta_PF_1/d_eta_PF` — arrays synced conditionally but no paired conditional `enter data`. If the `if` branch at sync time doesn't match the `if` at allocate time, this is undefined behaviour on the device.

---

### `src/core/MOM_continuity_PPM.F90`

No public derived types. Heavy kernel module: 80+ `do concurrent` loops across 3,123 lines, 44 target directives.

Local target-data regions (allocated/released within the same subroutine):
- `do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, uh_aux` — enter 1339, exit 1443. Mirror set for v at 1532/1614.
- Slope workspace `slp` — enter 2655/2798, exit 2740/2881.

**Sync points:** minimal. `Area_h` from/to at 254/278 (intermediate diagnostic reuse), `CS` update at 3067.

**Suspected issues:**
- Conditional `enter data` inside `if (CS%...)` branches at 1339, 1532, 1705. Safe only if the corresponding `do concurrent` is guarded by the same condition. Verify control flow before porting similar-shaped code.
- OBC segment arrays accessed via `OBC%segment(n)%...` inside concurrent loops at 1232–1246, 2215–2229. `OBC` is [CPU only] upstream — these rely on segment metadata being small and pre-synced by the caller. Audit.

---

### `src/core/MOM_CoriolisAdv.F90`

**Type: `CoriolisAdv_CS`** (lines 32–93)

| Field | Residency | Notes |
|---|---|---|
| `Coriolis_Scheme, KE_Scheme, PV_Adv_Scheme`, scalar configs | [GPU transient] | `update to(CS)` at 235 |
| `Time, diag` | [CPU only] | pointers |
| `id_*` (12 diag IDs) | [GPU transient] | scalars; synced with CS |

**Per-call target-data region (init at 247–316, teardown at 981–1005):** allocates ~20 workspace arrays — `Area_h, Area_q, dvdx, dudy, hArea_u/v, rel_vort, abs_vort, q, Ih_q, a, b, c, d, ep_u/v, KE, KEx, KEy`, plus conditional `dvSdx, duSdy, stk_vort, qS, uh_center, vh_center` under `Stokes_VF` or `Coriolis_En_Dis`. Correct lifecycle.

** Confirmed bugs in exit region (MOM_CoriolisAdv.F90:999–1005):**

```fortran
999: !$omp target exit data map(from: RV) if (CS%id_RV > 0)
1000: !$omp target exit data map(from: PV) if (CS%id_RV > 0) ! should be id_PV
...
1003: !$omp target exit data map(from: AD%rv_x_u) if (associated(AD%rv_x_u))
1004: !$omp target exit data map(from: AD%rv_x_u) if (associated(AD%rv_x_u)) ! duplicate
```

- Line 1000: guard uses `CS%id_RV` instead of `CS%id_PV`. `PV` copy-out is coupled to the wrong diag ID — if RV is requested but PV isn't, PV is still copied (benign); if PV is requested alone, it won't be (**functional bug — PV diag silently broken**).
- Line 1004: duplicates line 1003. The paired `AD%rv_x_v` exit is missing — that array never gets copied back, so any diag relying on `rv_x_v` is stale.

**Also suspicious:** line 316 `map(to: pbv, pbv%por_face_areaU, pbv%por_face_areaV) if (CS%Coriolis_En_Dis)` — if the flag is false the struct is not mapped, but downstream code may still reference `pbv`. Verify all `pbv` accesses in this file are guarded.

---

### `src/core/MOM_PressureForce*.F90`

`MOM_PressureForce.F90` is a thin dispatch layer — 2 directives (line 111 alloc FV CS, line 118 `update to(CS)`).

`MOM_PressureForce_FV.F90`: `PressureForce_FV_CS` mostly scalars/pointers; all field-level GPU state is transient workspace inside `PressureForce_FV_Boussinesq / _nonBoussinesq`. Nested target-data regions at 1267/1297–1350, 1425–1596, 1834–1950 allocate `pa, Z_0p, dpa, intx_dpa, inty_dpa, intz_dpa, EOSdom, press, T_int, S_int, rho_in_situ, dM, p0`. Conditional `Z_0p` at 1283 — if `use_EOS` is false, downstream reads at 1297 access unallocated device memory. ** Confirm `use_EOS` is always true in GPU builds, or gate the reads too.**

`MOM_PressureForce_Montgomery.F90`:

** Confirmed bug.** `PressureForce_Mont_CS` (lines 47–50) declares:

```fortran
real, allocatable :: PFu_bc(:,:,:)
real, allocatable :: PFv_bc(:,:,:)
```

Neither array appears in any `!$omp target enter data` across the file. They are the *outputs* of the Montgomery pressure-force computation. If Montgomery is selected and any part of its compute runs inside a target region writing to these arrays, it writes to unmapped memory.

**Current status:** Montgomery may not yet be on the GPU path in the current build — check with the user. But the fields are declared as if they were, which is a trap for future work.

Also: `SAL_CSp, tides_CSp` (55–56) are [CPU only]; if any Montgomery GPU kernel accessed them, it would segfault.

---

### `src/core/MOM_interface_heights.F90`, `src/core/MOM_variables.F90`

Zero target directives. Type definitions and CPU-only utilities. However, `MOM_variables.F90` defines the types (`thermo_var_ptrs, vertvisc_type, BT_cont_type, porous_barrier_type, accel_diag_ptrs, cont_diag_ptrs, ocean_internal_state, surface`) that are mapped *by other modules*. When porting consumers of these types, you inherit those modules' mapping discipline — check `MOM.F90` init and the owning module.

`MOM_interface_heights.F90` has 3 `do concurrent` loops but no target wrappers — those loops run on the host today. Candidate for porting if they sit on the hot path (they're called from pressure-force and dyn_split routines).

---

### `src/core/MOM.F90`

Top-level init/finalize — does most of the persistent mapping. See the "Big picture" box above.

**Issues:**
- Line 3016: `!!!$omp target enter data map(to: GV, GV%Rlay, GV%g_prime)` — commented out; enabled again at 3594. Indicates lifecycle uncertainty; worth a comment explaining why.
- Lines 2997–3008: maps `G` *and* its members separately. OpenMP struct mapping typically does not follow allocatables inside the struct, so explicit member maps are correct; but the `map(to: G)` at 2997 adds little beyond the scalar members. Consolidate or document.
- Line 672: `enter data map(to: forces, forces%taux, forces%tauy, forces%ustar)` with no matching `exit data` anywhere. ** Leaks device memory per run** (benign for a single run, problematic for in-process repeated inits e.g. coupled runs).
- Commented-out syncs at 853, 917–919, 1047–1049, 1110 — likely dead code from an earlier port iteration. Delete if truly unused (bit-rot risk).

---

### `src/parameterizations/lateral/MOM_hor_visc.F90`

**Type: `hor_visc_CS`** — ~40 fields, most `[GPU persistent]` via a large init block (lines 2894–3482). Full table in the agent report if needed; the noteworthy subset:

| Field group | Residency | Comment |
|---|---|---|
| `Kh_bg_xx/xy, Ah_bg_xx/xy, reduction_xx/xy, Laplac2_const_*, Biharm_const_*, Idx*, Idy*, dx2*, dy2*` | [GPU persistent] | Fully mapped |
| `Kh_Max_*, Ah_Max_*, Ah_Max_*_KS` | [GPU persistent, conditional] | Mapped under relevant scheme flag |
| **`n1n2_h, n1n1_m_n2n2_h, n1n2_q, n1n1_m_n2n2_q`** | [unmapped] | Anisotropic viscosity arrays. Allocated 2919–2922, **used in GPU regions at 1331–1338, never explicitly mapped.** If `anisotropic=true`, behavior is compiler-dependent. |
| **`m_const_leithy, m_leithy_max, Biharm6_const_xx/xy, Re_Ah_const_xx/xy`** | [unmapped] | Leith / Leith+E / Reynolds-bound arrays. Same issue — allocated conditionally, used in GPU regions, not mapped. Line 1377 even contains a `! TODO: !$omp target update from(...?)` comment. **Leith+E on GPU is incomplete.** |
| `Kh_bg_2d` | [CPU only] | Read from file at 2947, consumed at init in CPU loops at 3067, 3092–3099 to populate `Kh_bg_xx/xy` *before* those are mapped to device. Correct as-is. (Earlier agent report flagged this as a bug — verified and retracted.) |
| `diag` | [CPU only] | pointer |

**`do concurrent` kernels:** 56. **Sync points:** ~80 — many are genuine OBC roundtrips (lines 785–906), but a large fraction are ️ wasted roundtrips of the pattern:

```fortran
!$omp target update from(Kh) ! 1203
do j=...; Kh(i,j) = VarMix%Res_fn_h(i,j) * Kh(i,j); enddo
!$omp target update to(Kh) ! 1207
```

This pattern repeats for `Kh` at 1203–1207, 1212–1216, 1225–1248, 1284–1288, 1292–1299, and similarly for `Ah` at 1365–1500. Estimated 10–15 eliminable roundtrips — each one is a candidate for a 2-line replacement with `do concurrent`.

---

### `src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90`

**Type: `VarMix_CS`.** Arrays (`SN_u/v, L2u/v, cg1, Res_fn_*, Depth_fn_u/v, beta_dx2_*, f2_dx2_*, Rd_dx_h, slope_x/y, ebt_struct, sqg_struct, BS_struct, khth_struct, khtr_struct, kdgl90_struct, Laplac3_const_*, KH_u/v_QG`) are all declared `allocatable` but **no module-owned `enter data` directives exist** in this file. Only 2 syncs: `update from(h)` at 265 and 292 (before CPU `wave_speed` call).

**Question for the porting team:** is `VarMix` mapped by a *caller* at outer scope, or are `VarMix%Res_fn_h` etc. implicitly mapped when passed as a `type` to another module's target region? If the latter, it's working by accident. Making the lifecycle explicit here would prevent surprise breakage. `wave_speed_CS` (embedded at 190) is called CPU-side; keep [CPU only].

0 `do concurrent` kernels — this is a configuration/management module.

---

### `src/parameterizations/vertical/MOM_vert_friction.F90`

**Type: `vertvisc_CS`.**

| Field | Shape | Residency | Notes |
|---|---|---|---|
| `a_u, a_u_gl90, a_v, a_v_gl90` | 3D (nk+1) | [GPU persistent] | drag coefficients |
| `h_u, h_v` | 3D | [GPU persistent] | effective layer thickness |
| `a1_shelf_u, a1_shelf_v` | 2D pointers | [CPU only] | optional ice-shelf coupling |
| `kappa_gl90_2d` | 2D | [GPU transient] | optional GL90 (172) |
| Diagnostic arrays | various | [CPU only] | |

**`do concurrent` kernels:** 27. **Sync points:** 30+.

**Notable:** lines 674–677 map `ADp` struct plus its nested allocatables:

```fortran
!$omp target enter data map(to: ADp)
!$omp target enter data map(alloc: ADp%dv_dt_str, ADp%du_dt_str)
```

With some compilers, `map(to: ADp)` will create a device struct with dangling descriptors for allocatable members, then the second directive either refreshes them or conflicts. On nvfortran this is known-fragile; on ifx it's generally fine. **Worth a targeted test.**

Line 1436–1449 allocate `Ustar_2d, z_i, z_i_gl90, dz_harm, hvel, dz_vel, ...` on device with no visible matching exit in the same subroutine — exits are at 2059–2066 (~600 lines later). This is symmetric, just hard to read; consider a structured `target data` region.

---

### `src/parameterizations/vertical/MOM_set_viscosity.F90`

**Type: `set_visc_CS`** — mostly scalars ([CPU only]), plus `tideamp, bbl_u, bbl_v` (2D, [GPU transient]).

** Critical issue at lines 545–560:**

```fortran
do concurrent (i=is:ie, do_i(i,j) .and. (OBC%segnum_u(I,j) /= 0))
 call cvmix_kpp_composite_Gshape(sigma, Gat1, Gsig, dGdsig)
enddo
```

`cvmix_kpp_composite_Gshape` is a CPU-only external routine. Fortran 2018 allows calls inside `do concurrent` only to pure procedures, and GPU offload further demands `!$omp declare target`. Neither holds. On nvfortran with `-stdpar=gpu` or `-mp=gpu` this may compile but will execute on the host (silently losing parallelism) or crash. **Action:** either replace with inline math (cvmix `G` function is analytic and short), or hoist the loop out of `do concurrent`.

Also: the `tv%eqn_of_state` argument passed at line 374 to `calculate_density` — this is the EOS runtime-polymorphism hazard (see §EOS).

---

### `src/parameterizations/vertical/MOM_diabatic_aux.F90`

**Type: `diabatic_aux_CS`.** Scalar config + `createdH, penSW_diag, penSWflux_diag, nonpenSW_diag` ([GPU transient, diagnostic]).

Small, clean target region at 545–615: alloc workspace → one `teams loop` kernel → release. Only 3 sync points. Low risk.

`make_frazil` and `differential_diffuse_T_S` have no target directives — frazil is a candidate (pure arithmetic), the tridiagonal solver is genuinely hard to parallelize vertically.

---

### `src/tracer/MOM_tracer_advect.F90`

**Type: `tracer_advect_CS`.** Scalar-only: `dt, debug, useHuynhStencilBug, default_advect_scheme`, plus `diag` ([CPU only]) and `pass_uhr_vhr_t_hprev` (group_pass_type, [CPU only]).

**This file is heavily ported** — 70+ `do concurrent` kernels in 4 target regions (134–422, 512–824, 904–1250, 1298–1310). But it contains the most fragile pattern in the codebase:

 **Polymorphic registry traversal at lines 240, 400:**

```fortran
do m = 1, ntr
 !$omp target enter data map(to: Reg%Tr(m)%t)
enddo
```

`Reg` is a `tracer_registry_type, pointer`; `Reg%Tr` is an allocatable array of `tracer_type`; each `%t` is a 3D pointer into a tracer concentration field. The mapping relies on the compiler to:
1. Walk the pointer `Reg`
2. Index into the allocatable array `Tr(m)` on the host descriptor
3. Find and copy the 3D array `%t` at the device

Most OpenMP implementations do this incorrectly for deeply nested types. nvfortran in particular has had known bugs here. **Safer pattern:** flatten the registry once at init into a plain 4D array `tracer_conc(ntr, i, j, k)` (or pointer array of `device_ptr` with explicit `is_device_ptr` clauses), sync against it, and re-inject on exit.

 **Line 136:** `map(to: OBC)` for a structure with nested pointer arrays (`OBC%segment(n)%tr_Reg`). Same issue — what the OpenMP spec says will work is not what current compilers actually do reliably. Audit before trusting.

 **Line 385:** `update from(domore_k)` → CPU checks an iteration flag → re-enters GPU. Breaks device persistence and turns the iteration loop into ping-pong. If the iteration count is bounded, convert to a fixed-iteration loop with early-exit on device instead.

**Line 169:** `update from(local_advect_scheme)` — tiny, cheap; moves an integer back so the CPU can pick a stencil. Low-priority.

---

### `src/core/MOM_forcing_type.F90`

Big module (4462 lines), minimal GPU touch. Two types: `forcing` (thermodynamic) and `mech_forcing`, both ~50 `real, pointer, dimension(:,:)` surface-flux fields.

Only one target directive: `update to(U_star)` at line 1233, synced after CPU `find_ustar_fluxes`. The `do concurrent` loops at 1270–1295 are CPU-side — they are `do concurrent` for vectorization, not GPU offload (no surrounding target region).

**Implication for porting:** surface flux fields are [CPU only] today. Any GPU code that wants `forces%taux, forces%tauy, forces%ustar` needs an explicit `update to` before the region — as `MOM.F90:672` already does at init. Flux time-evolution over a step is CPU-side; if you port routines that *consume* fluxes inside a target region, the sync is still required.

`coupler_2d_bc_type` nested in `forcing%tr_fluxes` (line 249) is opaque w.r.t. GPU mapping — treat as [CPU only].

---

### EOS modules (`src/equation_of_state/`)

Files: `MOM_EOS.F90` (dispatcher), `MOM_EOS_base_type.F90` (abstract), and 9 concrete implementations (`_Wright, _Wright_full, _Wright_red, _linear, _Jackett06, _UNESCO, _Roquet_rho, _Roquet_SpV, _TEOS10`).

**The polymorphism problem, concretely.** The dispatcher in `MOM_EOS.F90` dispatches via type-bound procedures:

```fortran
call EOS%type%calculate_density_array(T, S, pressure, rho, is, npts, rho_ref=rho_ref)
! ^^^^^^^^^ polymorphic object (class(EOS_base))
```

`EOS%type` is `class(EOS_base), allocatable` — determined at init (`allocate(linear_EOS :: EOS%type)` or `allocate(buggy_Wright_EOS :: EOS%type)`, etc.). At runtime the call dereferences the v-table for the dynamic type and does an indirect branch.

**Why this is a blocker for GPU offload:**
- Indirect calls through a v-table require the v-table to be mapped to the device. For `class(EOS_base)` with 10 subclasses, each v-table entry points into host code. GPU runtimes do not reliably mirror v-tables; nvfortran in particular does not support this.
- Even if the vtable were mirrored, the GPU would need all 10 implementations compiled as device code and selected at runtime per element — costly and not what any current compiler actually does.

**Evidence this is already biting the port:** `MOM_EOS_Wright.F90` lines 106–113 contain this workaround:

```fortran
real elemental function density_elem_buggy_Wright(this, T, S, pressure)
 class(buggy_Wright_EOS), intent(in) :: this
 ! Body immediately calls a standalone function — `this` is discarded:
 density_elem_buggy_Wright = density_elem_buggy_Wright_loc(T, S, pressure)
end function
```

Comment at 104–105 explains: this is a known nvfortran bug where `class(...)` inside a GPU kernel produces runtime errors. The fix is manual inlining via a standalone function. **Any subclass that needs per-instance state (lookup tables, calibration) would silently break under this workaround** — the `this` is discarded.

**Recommended refactor.** Don't subclass at runtime. Replace with static dispatch resolved once at init:

```fortran
! At init (CPU):
select case (EOS%form_of_EOS)
 case (EOS_LINEAR) ; CS%eos_kind = EOS_LINEAR
 case (EOS_WRIGHT_FULL) ; CS%eos_kind = EOS_WRIGHT_FULL
 ...
end select

! At every call site (including inside `do concurrent`):
select case (CS%eos_kind)
 case (EOS_LINEAR) ; rho = density_elem_linear(T, S, p)
 case (EOS_WRIGHT_FULL) ; rho = density_elem_Wright_full(T, S, p)
 ...
end select
```

The `select case` is cheap (branch predictable, one scalar integer) and is fully static — the compiler inlines and GPU-compiles each branch. All EOS implementations already have standalone functions (the `_loc` variant in `MOM_EOS_Wright.F90` shows the pattern); the refactor is mostly renaming and adding `!$omp declare target` on the elemental functions.

**Alternative if refactor is deferred:** pre-compute density on the CPU once per step and pass the 3D field as input to GPU kernels. This works for code paths that only *read* density, but anything that computes density at many points per kernel (neutral diffusion, set-viscosity BBL, pressure force) pays a huge cost.

**EOS call sites inside / adjacent to GPU regions** (locations to audit):
- `MOM_PressureForce_FV.F90:1834–1950` — density in nested target region (currently likely compiled but relying on the `this`-discarding workaround).
- `MOM_PressureForce_Montgomery.F90:685–763` — same.
- `MOM_set_viscosity.F90:374, 732–738` — calculate_density in/near GPU region.
- `MOM_density_integrals.F90` — not yet GPU, but future target.
- `MOM_diabatic_driver.F90` — future target.
- `MOM_neutral_diffusion.F90` — future target; neutral-diffusion inner loop needs EOS at many points.

---

## Suspected bugs, prioritized

### Confirmed (file:line verified)

1. **`MOM_CoriolisAdv.F90:1000`** — `map(from: PV) if (CS%id_RV > 0)` — wrong guard; should be `CS%id_PV`.
2. **`MOM_CoriolisAdv.F90:1004`** — duplicate `map(from: AD%rv_x_u)`; line should be `AD%rv_x_v`. One-char typo. Silent diagnostic breakage.
3. **`MOM_PressureForce_Montgomery.F90:47–50`** — `PFu_bc, PFv_bc` allocatable, never mapped. If Montgomery PF runs on GPU, writes to unallocated device memory.
4. ~~`MOM_hor_visc.F90:2947` — `Kh_bg_2d` never mapped~~ **Retracted on verification.** The agent report was wrong. `Kh_bg_2d` is only consumed at init (lines 3067, 3092–3099) inside plain `do j / do i` loops *before* `CS%Kh_bg_xx/xy` are mapped to GPU. The values are baked into `Kh_bg_xx/xy`, which *are* mapped. `Kh_bg_2d` correctly stays CPU-only.
5. ~~`MOM_set_viscosity.F90:545–560` — `cvmix_kpp_composite_Gshape` in `do concurrent`~~ **Retracted on verification.** The function doesn't exist in this file. Lines 545–554 are plain array averaging (`T_vel, S_vel, Rml_vel, SpV_vel`) inside `do concurrent` — correct. The `do concurrent` at 560 with `DO_LOCALITY(local(k))` iterates OBC segments and contains only `if (OBC%segnum_u(I,j) > 0)` branches with plain assignments — also correct. Agent hallucination.

### Likely bugs (high confidence, lightly verified)

6. **`MOM_dynamics_split_RK2.F90:1601, 1608, 1614, 1657, 1669, 1677`** — `enter data` with no matching `exit data` in `end_dyn_split_RK2`. Orphans device memory for `taux_bot, tauy_bot, visc_rem_u/v, u/v_accel_bt, continuity_CSp, PressureForce_CSp, vertvisc_CSp, uhbt, vhbt`. One-run leak; multi-init coupled-model leak.
7. **`MOM.F90:672`** — `forces, forces%taux, forces%tauy, forces%ustar` mapped with no matching `exit data`. Same leak class.
8. **`MOM_barotropic.F90:1504, 1674–1676`** — conditional `update if(...) from(...)` for `uhbt0/vhbt0, eta_IC, eta_PF_1/d_eta_PF` where the `if` guard can diverge from the allocate-time guard. Latent; triggers only when the two predicates disagree.
9. **`MOM_tracer_advect.F90:240, 400`** — polymorphic registry mapping `Reg%Tr(m)%t` across nested pointer/allocatable descriptors. Compiler-dependent; flatten registry before GPU entry.
10. **`MOM_tracer_advect.F90:136`** — `map(to: OBC)` over a structure with nested pointer arrays (`segment(n)%tr_Reg`). Same class as #9.

### Probable gaps (need code inspection to confirm)

11. **`MOM_hor_visc.F90:1331–1338, 1377`** — anisotropic (`n1n2_*, n1n1_m_n2n2_*`) and Leith+E (`m_const_leithy, m_leithy_max, Biharm6_const_*, Re_Ah_const_*`) arrays used in GPU kernels but not mapped. File even has a `! TODO: !$omp target update from(...?)` at 1377. Leith+E on GPU is incomplete.
12. **`MOM_PressureForce_FV.F90:1283, 1297`** — `Z_0p` alloc conditional on `use_EOS`; downstream read at 1297 unconditional.
13. **`MOM_dynamics_split_RK2.F90`** — intent(inout) args (`h, u_inst, v_inst, uh, vh`) and grid masks (`G%mask2dCu/Cv`) used inside `do concurrent` without explicit `enter data` in this file. Works today because `MOM.F90` init covers them, but couples modules invisibly.
14. **`MOM_lateral_mixing_coeffs.F90`** — large transient workspace allocatables (`Res_fn_*, slope_x/y, ebt_struct, ...`) have no explicit mapping directives in this file. Rely on caller to map, but the contract is undocumented.
15. **`MOM_vert_friction.F90:674–677`** — `map(to: ADp)` followed by `map(alloc: ADp%du_dt_str, ADp%dv_dt_str)`. Order-sensitive and compiler-fragile for nested allocatables.

### ️ Wasted roundtrips (performance, per-step)

16. **`MOM_dynamics_split_RK2.F90`**: 597/599, 797/800, 848/858, 984/986, 1093/1096 — OBC zeroing + vertFPmix/vertvisc. Eliminable by porting the small CPU routine.
17. **`MOM_hor_visc.F90:1203–1260, 1365–1500`** — ~10–15 `Kh`/`Ah` from/CPU-multiply/to roundtrips. Each replaceable with a 3-line `do concurrent`.
18. **`MOM_dynamics_split_RK2.F90:767–768, 1070–1072`** — `update from(up, vp)` / `update from(u_inst, v_inst)` solely to zero a CPU-local array (`uold(:,:,:) = 0.0`). Gratuitous — no actual data is needed from the GPU. Delete the `update from`.

---

## Suggested porting order (practical)

If the goal is to increase on-device dwell time and reduce sync per step, tackle in this order:

1. **Fix the five confirmed bugs (§1–5)** — they are cheap and unblock trust in the port.
2. **Close the `exit data` gaps (§6–7)** — low-effort, removes a leak class, makes re-init coupled-mode safe.
3. **Kill the wasted roundtrips in `MOM_hor_visc` (§17)** — ~10–15 trivial `update from/to` pairs replaceable with `do concurrent`. Biggest per-step ratio win.
4. **Port `Stokes_PGF` and `open_boundary_zero_normal_flow`** — small CPU routines sitting on the hot path in `MOM_dynamics_split_RK2`. Each one killed removes 2–4 sync points per step.
5. **Port `vertFPmix`** — bigger job, but lines 797/800 and 1093/1096 force 5-array sync pairs per step. `vertvisc` is much bigger; `vertFPmix` alone pays off.
6. **EOS refactor (§EOS)** — blocks Montgomery, neutral diffusion, and anything inside set_viscosity's inner EOS calls. Best done before pushing further into set_viscosity and density integrals.
7. **Tracer registry flattening (§9)** — required before trusting `MOM_tracer_advect.F90` under broader compiler coverage. Can be deferred if nvfortran is the only target and current tests pass, but it's a ticking correctness bomb.
8. **hor_visc gaps (§11)** — only matters under anisotropic or Leith+E; re-prioritize if those schemes are used in production.

---

## Known-CPU subroutines that are called from GPU context (sync boundaries)

These are the walls you hit when porting more of the dynamics core:

- `Stokes_PGF` (MOM_wave_interface)
- `open_boundary_zero_normal_flow`, `radiation_open_bdry_conds`, `update_OBC_data` (MOM_open_boundary, MOM_boundary_update)
- `CorAdCalc` (MOM_CoriolisAdv has GPU regions inside, but the top-level call is CPU — the per-call transient arrays live on device, but the dispatch is host-side)
- `vertFPmix`, `vertvisc` (MOM_vert_friction — vert_friction has GPU regions but these entry points have not been fully migrated)
- `wave_speed` (MOM_lateral_mixing_coeffs callers)
- `thickness_to_dz` (MOM_interface_heights)
- All EOS dispatch until refactored (§EOS)
- `pass_var`, `pass_vector`, `do_group_pass`, halo passes (MPI — genuinely CPU, unavoidable without GPU-aware MPI)

---

---

# Unmapped array inventory (verified)

Dedicated audit pass to answer *"what CS-struct array fields are latent GPU requirements we haven't mapped yet?"* Three parallel audits were run; their raw outputs contained several false positives (details at the end of this section), so every claim below has been verified by direct `grep` against source.

## Confirmed unmapped fields used in GPU regions (new bugs)

All of these are `allocatable`/`pointer` array fields declared in a CS struct, used inside a `!$omp target` region or `do concurrent`, with **no matching `enter data`** anywhere in their owning file. If the code path that uses them fires, the device reads uninitialized memory.

### `MOM_hor_visc.F90` — `hor_visc_CS`

The existing §map already flagged this as "probable gaps"; the audit confirms with line numbers. These are the unmapped 2D CS fields consumed by GPU kernels:

| Field | Decl line | Consumer (in target region) | Trigger condition |
|---|---|---|---|
| `Kh_Max_xx` | 155 | 1214, 1267, 1277 | `CS%better_bound_Kh` or anisotropic bound |
| `Ah_Max_xx` | 156 | 1346, 1434, 1460, 1464 | `CS%better_bound_Ah` |
| `Ah_Max_xx_KS` | 157 | 1472, 1474 | `use_Kh_bg_2d`-adjacent biharmonic bound |
| `n1n2_h` | 159 | 1255, 1336 | `CS%anisotropic` |
| `n1n1_m_n2n2_h` | 160 | 1336 | `CS%anisotropic` |
| `Kh_Max_xy`, `Ah_Max_xy`, `n1n2_q`, `n1n1_m_n2n2_q` | 175–180 | q-point analogs of above | same as above |
| `m_const_leithy` | 187 | 1388 | `CS%use_Leithy` |
| `m_leithy_max` | 188 | 1391 | `CS%use_Leithy` |
| `Biharm6_const_xx`, `Biharm6_const_xy` | 204, 211 | 1370, 1384, 1407 (+ q analogs) | `CS%Leith_Ah` with biharmonic |
| `Re_Ah_const_xx`, `Re_Ah_const_xy` | 208, 215 | (Reynolds-number bounded Ah) | `CS%Re_Ah > 0` |

**Fix recipe** for each group — conditional enter at the end of `hor_visc_init`, mirrored exit in `hor_visc_end`:
```fortran
if (CS%anisotropic) then
 !$omp target enter data map(to: CS%n1n2_h, CS%n1n1_m_n2n2_h, CS%n1n2_q, CS%n1n1_m_n2n2_q)
endif
if (CS%use_Leithy) then
 !$omp target enter data map(to: CS%m_const_leithy, CS%m_leithy_max)
endif
if (CS%Leith_Ah) then
 !$omp target enter data map(to: CS%Biharm6_const_xx, CS%Biharm6_const_xy)
endif
! ... matching Kh_Max_*/Ah_Max_*/Re_Ah_* per their own guards
```

### `MOM_variables.F90` — `porous_barrier_type`

| Field | Decl line | Consumer | Status |
|---|---|---|---|
| `por_face_areaU` | 358 | `MOM_CoriolisAdv.F90:316` → `enter data map(to: pbv%por_face_areaU) if (CS%Coriolis_En_Dis)` | conditionally mapped |
| `por_face_areaV` | 359 | same | conditionally mapped |
| **`por_layer_widthU`** | 360 | `MOM_set_viscosity.F90:1001` (inside an `if (use_BBL_EOS)` branch) | ** unmapped** |
| **`por_layer_widthV`** | 361 | `MOM_set_viscosity.F90:1002` | ** unmapped** |

`MOM_porous_barriers.F90` writes these fields CPU-side every step (lines 241–271), after which they'd need to be synced up to device. When set_viscosity is ported beyond its current stub, `por_layer_widthU/V` need `enter data` mirroring what `por_face_areaU/V` already do, plus a per-step `update to` after the porous barrier recompute.

### `MOM_PressureForce_Montgomery.F90` — `PFu_bc / PFv_bc`

Already in the §"Confirmed bugs" list above (item #3). Re-flagged here because it fits the exact same pattern: allocatable array field, no `enter data`, consumed by an eventual GPU kernel.

## False-positive retractions (audit hallucinations)

To prevent them from re-appearing in future audits:

- **`barotropic_CS%frhatu/frhatv`** — claimed unmapped; mapped at `MOM_barotropic.F90:6238`. Explicit `update from` also present at line 1738 (post-diag copy). Mapped.
- **`barotropic_CS%eta_cor`** — claimed unmapped; mapped at `MOM_barotropic.F90:6239`. Written CPU-side (`set_dtbt` at 5411, 5416, filled at init 5879), then `update to` applies at 6239. Mapped.
- **`MOM_dyn_split_RK2_CS%PFu_Stokes, PFv_Stokes`** — claimed used inside target region; **not inside a target region**. Line 537 is a CPU call to `Stokes_PGF`; lines 543–552 are plain `do` loops (not `do concurrent`). Correctly CPU-only. They sit between an `update from` (line 535) and a later `update to` that resyncs the CPU-side modified `CS%PFu/PFv`. No bug.
- **`vertvisc_type%Ray_u/Ray_v`** — claimed unmapped; conditionally mapped at `MOM_vert_friction.F90:740, 939`, released at 1106–1107.
- **`vertvisc_type%Kv_shear, Kv_slow, TKE_turb, Kd_shear, bbl_thick_*, kv_bbl_*, ustar_BBL`** — claimed unmapped GPU-consumers. **Verified: only consumer in GPU regions is `post_data` at `MOM_vert_friction.F90:2080–2081`, which is CPU-side.** These fields are populated by `MOM_kappa_shear`, `MOM_set_diffusivity`, `MOM_set_viscosity` — none of which are yet GPU-ported. So they're correctly [CPU only] *today*; they will become [GPU persistent] when those modules are ported (Tier 2). Not a current bug; add to pre-flight list.

## `MOM_remapping` — second `class(...)` polymorphism blocker

**This was not in the existing §EOS section and deserves equal weight.** Verified at `src/ALE/MOM_remapping.F90:82`:

```fortran
class(Recon1d), pointer :: reconstruction => Null()
```

`src/ALE/Recon1d_type.F90` defines the abstract base; `src/ALE/Recon1d_*.F90` has **17 concrete subclasses**: `PCM, PLM_CW, PLM_CWK, PLM_hybgen, PPM_CW, PPM_CWK, PPM_H4_2018, PPM_H4_2019, PPM_hybgen, MPLM_CWK, MPLM_WA, MPLM_WA_poly, EPPM_CWK, EMPLM_CWK, EMPLM_WA, EMPLM_WA_poly` (plus the abstract type). Dispatch is via `select type` in the remapping hot loop.

**Implication for the Tier 1 ALE port:** porting `MOM_remapping.F90` (issue #119) and `MOM_ALE.F90` (issue #120) hits exactly the same GPU wall as EOS. Needs a parallel static-dispatch refactor: introduce `CS%recon_kind` integer (matching the existing `REMAPPING_PLM_CW`, `REMAPPING_PPM_H4_*` enum values likely already present in `regrid_consts`), replace every `select type(this%reconstruction)` with `select case(CS%recon_kind)` calling standalone module-level kernels. Each Recon1d_*.F90 already has the numeric kernels written as methods; the refactor is mostly mechanical renaming + adding `!$omp declare target`.

**Update to Tier 0:** the Tier 0 refactor isn't just "EOS static dispatch" — it's "runtime-polymorphism-to-static-dispatch, applied to both EOS and Recon1d." They're independent code paths (EOS in `src/equation_of_state/`, Recon1d in `src/ALE/`), can be done in parallel PRs. Both block Tier 1.

## Pre-flight inventory for Tier 1/2 modules

Summaries of CS array fields that will need `enter data` when these modules are ported. Array-size notes use 360×180×75 (roughly 1° × 1° ACCESS-OM3 shape); 8-byte reals. The earlier audit's GB estimates were unit-confused — corrected below.

### `MOM_neutral_diffusion.F90` — `neutral_diffusion_CS`

Largest field inventory of any unported module. All `[persistent]`:

| Field | Dims | Bytes (per field) |
|---|---|---|
| `hbl`, `ns` | 2D (i,j) | ~520 KB |
| `Coef_h` | 3D (i,j,k+1) | ~40 MB |
| `uPoL, uPoR, uHeff` | 3D (IB,j,k+1) | ~40 MB each |
| `vPoL, vPoR, vHeff` | 3D (i,JB,k+1) | ~40 MB each |
| `uKoL, uKoR, vKoL, vKoR` | 3D integer | ~20 MB each |
| `dRdT, dRdS, Tint, Sint, Pint` | 3D | ~40 MB each |
| `ppoly_coeffs_T, ppoly_coeffs_S` | 4D (i,j,k,order≈3) | ~120 MB each |
| `T_i, S_i, P_i, dRdT_i, dRdS_i` | 4D (i,j,k+1,order) | ~120 MB each |
| `stable_cell` | 3D logical | ~5 MB |

**Total persistent device footprint: ~1.5 GB** (not the 3–4 TB the audit claimed — arithmetic error, `180*90*76*8 ≈ 10 MB` not `10 GB`). Comfortable on any modern GPU.

Embeds `remap_CS` (line 120) — brings the `Recon1d` polymorphism issue above into scope for neutral_diffusion too. External pointers `EOS, KPP_CSp, energetic_PBL_CSp` stay CPU or pointer-flat.

### `MOM_thickness_diffuse.F90` — `thickness_diffuse_CS`

| Field | Dims | Size |
|---|---|---|
| `GMwork` | 2D | ~520 KB |
| `Kh_eta_u, Kh_eta_v` | 2D (face) | ~520 KB each |
| `KH_u_GME, KH_v_GME` | 3D (face, k) | ~40 MB each |
| `khth2d` | 2D (per-step) | ~520 KB |
| `diagSlopeX, diagSlopeY` | 3D (per-step diagnostic) | ~40 MB each |

**Total: ~160 MB persistent + 80 MB diagnostic transient.** Unlikely to be a constraint.

### `MOM_energetic_PBL.F90` — `energetic_PBL_CS`

All modest: `ML_depth, BBL_depth` are 2D (~520 KB each), `shape_function` is 1D (~1 KB), `ML_c` is size-18 compile-time. Total &lt;2 MB persistent. Easy to map.

### `MOM_mixed_layer_restrat.F90` — `mixedlayer_restrat_CS`

All 2D: `MLD_filtered, MLD_filtered_slow, wpup_filtered, MLD_Tfilt_space, Cr_space`. ~520 KB each, ~3 MB total. Trivial to map.

### `MOM_kappa_shear.F90`, `MOM_isopycnal_slopes.F90`, `MOM_set_diffusivity.F90`, `MOM_density_integrals.F90`

**No persistent arrays in their CS structs** — all working memory is local to the subroutines (per-call transient). Mapping discipline will be at the routine level, not the CS level. `set_diffusivity_CS` carries pointers to *other* modules' CS structs (`Kappa_shear_CS`, `CVMix_*_cs`, `int_tide_CS`, embedded `tidal_mixing_cs`); port ordering matters because each pointed-to CS needs its own mapping work first.

### `MOM_regridding.F90` — `regridding_CS`

Only short 1D arrays (`coordinateResolution, target_density, max_interface_depths, max_layer_thickness` — all ~1–2 KB). Plus conditional pointer-to-subtype structs (`zlike_CS, sigma_CS, rho_CS, hycom_CS, adapt_CS, hybgen_CS`) — each grid-type backend has its own small workspace. Negligible size; the port challenge is the coordinate-specific dispatch (same pattern as EOS/Recon1d but only 5 subtypes).

### `MOM_ALE.F90` — `ALE_CS`

Embeds `regridCS, remapCS, vel_remapCS` — so it inherits the remapping CS field set above. `hybgen_unmixCS` is conditional. Diagnostic ID arrays (`id_tracer_remap_tendency, do_tendency_diag`) are 1D integers/logicals — tiny. No novel mapping surface beyond what `remap_CS` requires.

## Audit-quality note

Running three agents in parallel produced overlapping coverage but also three classes of error worth recording:

1. **False positives on well-mapped fields** (audit 1: `frhatu`, `frhatv`, `eta_cor`). Caused by agents searching only within a narrow line range rather than the whole file. Countermeasure: always verify an "unmapped" claim with an unconstrained grep for `enter data.*FIELD` before trusting.
2. **False positives on CPU-only-by-design fields** (audit 1: `PFu_Stokes`, `Ray_u/v`, `Kv_shear`). Caused by agents not distinguishing `do concurrent` (GPU-eligible) from plain `do` (CPU), and not checking whether consumers are actually yet-to-be-ported. Countermeasure: for any claimed bug, verify the consumer is inside a `!$omp target` scope today.
3. **Unit/arithmetic errors** (audit 3: reported `180×90×76×8 bytes = 87 GB`, actual ≈ 10 MB; over by 10000×). Countermeasure: sanity-check any memory-pressure claim against `n_i * n_j * n_k * 8` before citing.

Confirmed findings in this section have all been grep-verified at the file:line level.

---

*This map reflects the state of the `dev/gpu` branch at commit `6474597b1`. Regenerate when the branch has moved more than a few hundred lines in the dynamics/parameterization areas.*

---

# Porting plan for the rest of the dynamical core + ALE

Scope: everything in `issues.json` labelled `gpu` + `ACCESS-OM3 config`, restricted to dynamics, EOS, ALE/remap, and their dependencies. Sorted by **dependency tier**, then by **ACCESS-OM3 benchmark cost** (the seconds reported in each issue body — treat as rough indicator of hot path weight).

**Reading key per row:** `issue # — file / routine (cost) — blockers — approach — notes`.

The "cost" column is total ACCESS-OM3 time when available, otherwise benchmark time. It's a routing signal, not a target: the `mixedlayer_restrat_bodner` 365s and `int_density_dz_generic_plm` 977s numbers dominate everything else and should shape the order.

---

## Tier 0 — Unblock the EOS and its scalar math primitives

**Why first:** nearly every hot routine in Tier 1 calls `calculate_density`, `calculate_density_derivs`, or `calculate_spec_vol`. With runtime polymorphism in place, porting the callers produces either compile errors or silent host fallback. Fix this once, unblock ~20 downstream issues.

| # | File / routine | Cost (s) | Approach |
|---|---|---|---|
| #72, #88 | `MOM_EOS_base_type.F90`, `MOM_EOS.F90` dispatcher | (structural) | **Replace polymorphic calls with static dispatch** as described in §EOS. Introduce `CS%eos_kind` integer at init (already exists as `form_of_EOS`), replace every `EOS%type%calculate_*` call with a `select case` over `eos_kind` invoking module-level standalone functions (the `_loc` pattern already used in `MOM_EOS_Wright.F90:106–113`). Add `!$omp declare target` to each elemental and each standalone function. Audit gptl points: `calculate_density_{1d,derivs_1d,derivs_array}`, `calc_spec_vol_derivs_1d`, `EOS_domain`. |
| #115 | `MOM_EOS_Roquet_rho.F90` — `density_anomaly_elem`, `density_elem`, `calculate_density_derivs_elem`, `calculate_density_array`, `calculate_specvol_derivs_elem` | **6246 + 898 + 788 + 187 + 8** — the single biggest hot spot in ACCESS-OM3 | Per @edoyango, upstream PR #156 addresses this. Rebase/cherry-pick when the dispatch refactor lands; Roquet is ACCESS-OM3's default EOS so it dominates. These are pure math; once `!$omp declare target` is on the elementals and the dispatcher no longer needs `class(*)`, they should port cleanly. |
| #122 | `MOM_temperature_convert.F90` — `potemp_to_constemp`, `constemp_to_potemp`, `dtc_dtp`, `potemp_to_constemp` | 31.9 + 16.5 + 11.5 + 3.95 | Pure elemental math on scalars. Annotate with `!$omp declare target` and wrap call sites in `do concurrent`. No control-flow issues. Low effort. |
| #128 | `MOM_TFreeze.F90` — `calculate_tfreeze_teos_poly_array` | 18.4 | Same pattern as #122. Add `declare target`, make callable inside GPU regions. |
| #124 | `MOM_intrinsic_functions.F90` — `cuberoot`, `rescale_cbrt_`, `descale_`, `invcosh` | 35.9 + 4.34 + 3.47 | Pure math. `declare target` + inline into `do concurrent`. Note: `cuberoot` is probably called from deep inside EOS/Roquet; confirm it's on the GPU call path before bigger EOS work. |
| #72 | `EOS_base_type.F90` — `a_calculate_density_derivs_array`, `a_calculate_specvol_derivs_array` | 0.32 + 1.36 | Default implementations in the abstract base. Once dispatch is static, these default loops become regular code — annotate and run as `do concurrent` over `is:ie`. |

**Expected outcome:** after Tier 0, any `do concurrent` in a caller can call density/spec-vol as a regular function. Unblocks Tier 1 and Tier 2.

Reference: closed issue #31 in `issues.json` documents a prior "EOS data transfer error in Nvidia 25.5" — useful prior art when touching this area.

---

## Tier 1 — ALE / remap / regridding cluster (user has context here)

**Why group together:** these are called back-to-back in a step (`ALE_regrid` → `ALE_remap_*` → `remap_*` → `interpolate_column`). They share the same `remapping_CS` and work on the same per-column reconstruction data, so device residency should be designed across the whole cluster. Keep regrid-remap data structures on device end-to-end; the CPU only needs to orchestrate.

| # | File / routine | Cost (s) | Approach |
|---|---|---|---|
| #120 | `MOM_ALE.F90` — `ale_plm_edge_values`, `ale_remap_tracers`, `ale_remap_velocities`, `ale_remap_set_h_vel`, `ale_remap_scalar`, `ale_regrid` | 383 + 23 + 15 + 8.7 + 2.5 + 1.8 | User self-assigned. Strategy: (a) `enter data` the `ALE_CS` and its per-column working arrays at init; (b) port `ale_plm_edge_values` first — it's embarrassingly parallel in (i,j) with a k-column serial inner loop, maps cleanly to `do concurrent (j,i)` with inner `do k`; (c) `ale_remap_*` wraps `remapping_core_h` — needs #119 done in sync. |
| #119 | `MOM_remapping.F90` — `remap_src_to_sub_grid`, `average_value_ppoly`, `intersect_src_tgt_grids`, `remap_sub_to_tgt_grid_om4`, `reintegrate_column`, `interpolate_column`, `build_reconstructions_1d`, `remapping_core_h` | 204 + 157 + 109 + 101 + 41 + 7.65 + 6.3 + 5.5 | User self-assigned `interpolate_column`. These are **column-local** kernels — each column is independent, so the outer loop maps to `do concurrent (j,i)` and inner column work runs serially on device. Need to pass `remapping_CS` into GPU regions: likely needs field-level `enter data` for its reconstruction workspaces. |
| #97 | `MOM_regridding.F90` — `get_rho_cs`, `build_zstar_grid`, `calc_h_new_by_dz`, `adjust_interface_motion`, `filtered_grid_motion`, `regridding_main`, `set_h_neglect`, `inflate_vanished_layers_old` | 12.2 + 10.5 + 8.0 + 4.6 + 1.5 + 1.3 + 0.42 + 0.20 | Called once per regrid step. `regridding_main` dispatches by grid kind — refactor to static dispatch similar to EOS (fewer subclasses, easier). z*, sigma, rho grids each need their inner kernel ported. |
| #75 | `MOM_isopycnal_slopes.F90` — `calc_isoneutral_slopes`, `vert_fill_ts` | 74.0 + 0.23 | `calc_isoneutral_slopes` calls `calculate_density_derivs` inside a 3D loop — **blocked on Tier 0**. Once EOS is static, this is a straightforward 3D `do concurrent`. |
| #117 | `MOM_neutral_diffusion.F90` — `neutral_diffusion`, `neutral_diffusion_calc_coeffs`, `neutral_surface_flux`, `ppm_edge_`, `find_neutral_surface_positions_continuous`, `ppm_ave`, `plm_diff`, `signum_`, `ppm_left_right_edge_values_`, `interface_scalar`, `fv_diff_` | 754 + 228 + 166 + 55 + 54 + 54 + 41 + 22 + 18 + 12 + ... | Hot. Neutral diffusion needs density repeatedly per cell-face — **the canonical case where EOS static dispatch pays off.** Strategy: port `neutral_diffusion_calc_coeffs` first (it builds the coefficients; everything else consumes them). The PPM helpers (`ppm_edge_`, `ppm_ave`, `plm_diff`) are pure, short — annotate with `declare target`. |

**Shared concern:** `remapping_CS`, `regridding_CS`, `neutral_diffusion_CS` — each should be inspected for allocatable working arrays that need explicit `map(alloc:)` at init. If any contains polymorphic or pointer-to-pointer fields, they need flattening (same class as the tracer registry issue in §MOM_tracer_advect.F90).

---

## Tier 2 — Diabatic column physics + mixing

**Why grouped:** all called from `diabatic_driver` and `diabatic_ale` (23.7s) on the column axis. Most have a (j,i) outer loop over vertical-column kernels. EPBL and KPP column solvers contain serial k-loops with intra-column dependencies.

| # | File / routine | Cost (s) | Approach |
|---|---|---|---|
| #66 | `MOM_diabatic_driver.F90` — `layered_diabatic`, `diabatic_ale`, `diabatic`, `extract_diabatic_member`, `diabatic_driver_init`, `diabatic_ale_legacy` | 23.7 + 0.73 + 0.13 | Orchestrator; mostly calls sub-routines. Port piecewise, following the order below (set_viscosity → set_diffusivity → EPBL → kappa_shear → opacity → applyboundaryfluxesinout → tracer_vertdiff → regrid/remap). |
| #121 | `MOM_energetic_PBL.F90` — `epbl_column`, `find_pe_chg_orig`, `energetic_pbl`, `find_pe_chg_orig_`, `mstar_langmuir`, `find_mstar` | 91.6 + 31.6 + 16.3 + 9.5 + 0.5 + 0.5 | `epbl_column` is a classical column solver: parallelize (j,i), keep k serial. Many iterative scalar solves inside — check if any call `cuberoot` or EOS (if so, Tier 0 blocks). |
| #59 | `MOM_kappa_shear.F90` — `calc_kappa_shear_vertex`, `find_kappa_tke`, `calculate_kappa_shear`, `kappa_shear_column`, `calculate_projected_state`, `kappa_shear_init`, `kappa_shear_is_used` | 82.0 + 1.23 + 0.41 + 0.23 + 0.12 | Per @edoyango's comment in #59: needs refactor. The "entry point" is parallel over (i,j) but `kappa_shear_column` has a serialised k-loop with iterative convergence (`do while`). Options: (a) accept the serialized k-loop on device (fine if work per column is substantial), (b) refactor to fixed max-iterations, kill ping-pong updates. EOS-dependent — Tier 0 blocks. |
| #68 | `MOM_set_viscosity.F90` — `find_l_open_convex`, `find_l_open_concave_trigonometric`, `set_viscous_ml`, `set_viscous_bbl`, `set_u_at_v`, `set_v_at_u` | 64.3 + 8.7 + 0.26 + 0.22 + ... | Per @edoyango: in progress. Known blocker: EOS dependency in `set_viscous_ml` (Tier 0). `find_l_open_convex/concave` are the hottest — pure arithmetic root-finders over (i,j), should port cleanly. The `do concurrent` blocks at 545–554 already handle T/S/SpV averaging and OBC-segment iteration; the harder loops are the geometry root-finders. |
| #70 | `MOM_set_diffusivity.F90` — `find_n2`, `set_diffusivity`, `find_tke_to_kd`, `set_bbl_tke`, `add_drag_diffusivity`, `set_density_ratios`, `add_lotw_bbl_diffusivity` | 1.27 (lotw) + 0.20 (n2) + 0.16 (bbl) + 0.12 + ... | Issue body flags `set_diffusivity` needs refactor: *"jki loop with routine calls"* — the routine calls inside that order probably aren't `declare target`; inline/copy out. `set_density_ratios` is EOS-dependent (Tier 0). |
| #87 | `MOM_bkgnd_mixing.F90` — `calculate_bkgnd_mixing` | 0.04 | Tiny. Just annotate and wrap in `do concurrent`. |
| #126 | `MOM_tidal_mixing.F90` — `add_int_tide_diffusivity`, `setup_tidal_diagnostics`, `tidal_mixing_h_amp` | 22.1 + 2.5 + 0.12 | Column-local. Port after #70 — same structure. |
| #77 | `MOM_opacity.F90` — `opacity_from_chl`, `absorbremainingsw`, `extract_optics_slice`, `set_opacity` | 17.2 + 0.12 + 0.05 + 0.03 | Column-local; no EOS. Straightforward 3D `do concurrent`. |
| #123 | `MOM_full_convection.F90` — `full_convection`, `is_unstable_`, `smoothed_drdt_drds` | 31.4 + 11.0 + 7.06 | Column k-loop with inter-level dependency. EOS calls (density derivs) inside — Tier 0 blocks. |
| #69 | `MOM_diabatic_aux.F90` — `applyboundaryfluxesinout`, `find_uv_at_h`, `adjust_salt`, `set_pen_shortwave`, `make_frazil`, `tridiagts_eulerian` | 25.8 + 1.11 (ALE) + 0.6 + 0.3 + 0.11 + 0.15 (ALE) | From §map: file has clean target-data region already (lines 545–615). Extend pattern to `applyboundaryfluxesinout` (hot), `find_uv_at_h` (trivial 3D), `make_frazil` (k-column serial). `tridiagts_eulerian` has a tridiagonal solve along k — serialize k, parallelize (j,i). |
| #83 | `MOM_tracer_diabatic.F90` — `applytracerboundaryfluxesinout`, `tracer_vertdiff_eulerian`, `tracer_vertdiff` | 12.2 + 5.8 + 0.09 | Tridiagonal per column. Same (j,i) parallel, k-serial pattern. |
| #142 | `MOM_diagnose_MLD` — `diagnosemldbydensitydifference` | 1.34 | EOS call inside loop (issue body flags). Blocked by Tier 0. |

---

## Tier 3 — Lateral transport + mixing coefficients

Most of these are already partially ported — see existing §MOM_hor_visc and §MOM_lateral_mixing_coeffs sections for known gaps. The issues below add concrete hot paths not yet covered.

| # | File / routine | Cost (s) | Approach |
|---|---|---|---|
| #81 | `MOM_mixed_layer_restrat.F90` — `mixedlayer_restrat_bodner`, `mu_`, `rmean2ts_`, `mixedlayer_restrat`, `mixedlayer_restrat_bml` | **365.5** + 44.7 + 1.14 + 0.16 + 0.13 | **Second-biggest hot spot overall.** Bodner scheme dominates ACCESS-OM3. Port whole file. Calls `mu_` and `rmean2ts_` heavily — those need `declare target`. Inspect for EOS dependency. |
| #62 | `MOM_thickness_diffuse.F90` — `streamfn_solver`, `thickness_diffuse_full`, `thickness_diffuse` | 25.4 + 1.39 + 0.12 | `thickness_diffuse_full` has density derivs inside — Tier 0 blocks. `streamfn_solver` is a vertical solve per (i,j) column — standard pattern. |
| #61 | `MOM_tracer_hor_diff.F90` — `tracer_epipycnal_ml_diff`, `tracer_hordiff`, `tracer_hor_diff_init` | 1.37 + 0.27 + 0.01 | `tracer_epipycnal_ml_diff` is EOS-dependent. `tracer_hordiff` is standard lateral diffusion, no EOS — port straight. |
| #74 | `MOM_lateral_mixing_coeffs.F90` — `calc_visbeck_coeffs_old`, `calc_slope_functions`, `calc_depth_function`, `calc_slope_functions_using_just_e`, `calc_resoln_function`, `varmix_init` | 33.7 + 0.05 + 0.04 + 0.23 + 0.015 + 0.005 | Per §map, this file has workspace allocatables without explicit mapping. Before porting, add the missing `enter data`/`exit data` for `Res_fn_*, slope_x/y, ebt_struct, ...`. `calc_visbeck_coeffs_old` is the hot one. |
| #71 | `MOM_wave_speed.F90` — `tdma6`, `wave_speed`, `tridiag_det` | 136.4 + 0.3 + 0.085 | `tdma6` is the hottest tridiagonal solver in the code. Serial k per column, parallel (j,i). EOS-dependent (Tier 0 blocks). |
| #129 | `MOM_MEKE.F90` — `step_forward_meke`, `meke_lengthscales_0d`, `meke_lengthscales`, `meke_equilibrium` | 11.1 + 3.0 + 1.31 + 0.38 | 2D fields, no EOS. Port straight. |
| #118 | `MOM_density_integrals.F90` — `int_density_dz_generic_plm`, `int_density_dz_generic_pcm`, `int_density_dz` | **977.8** + 5.9 + 0.11 | Third-biggest hot spot. Per @edoyango: `int_density_dz_generic_pcm` in progress for Roquet. Both PLM and PCM variants need static EOS dispatch. |
| #60 | `MOM_PressureForce_FV.F90` — `pressureforce_fv_bouss` | 1.94 | Issue says "mostly done, but requires refactor of EOS" — Tier 0 blocks final completion. |
| #73 | `MOM_PressureForce_Montgomery.F90` — `set_pbce_bouss` | 0.29 | EOS-dependent + the `PFu_bc/PFv_bc` never-mapped bug from §map. Fix the mapping first, then port once Tier 0 lands. |

---

## Tier 4 — Support / utilities (do when convenient)

| # | File / routine | Cost (s) | Notes |
|---|---|---|---|
| #63 | `MOM_tracer_advect.F90` — `advect_y`, `advect_x`, `advect_tracer` | 0.73 + 0.35 + 0.10 | Already heavily ported. Outstanding issue: polymorphic registry traversal at lines 240, 400 (see §map #9). **Don't port more callers until registry is flattened.** |
| #65 | `MOM_interface_heights.F90` — `thickness_to_dz_3d` (), `find_eta_3d`, `find_eta_2d`, `thickness_to_dz_jslice`, `dz_to_thickness_tv`, `find_rho_bottom` | 0.68 + 0.10 + 0.04 + 0.015 + 0.005 + 0 | `thickness_to_dz_3d` done (flag-gated). Issue body: callers of `_jslice` should switch to the 3D interface. `find_rho_bottom` has EOS call. |
| #78 | `MOM_forcing_type.F90` — `fluxes_accumulate`, `myalloc_2d`, `copy_back_forcing_fields`, `mech_forcing_...` | 11.5 + 0.69 + 0.42 | §map says this file is mostly CPU — that's correct and intentional. `fluxes_accumulate` is the hot one, but it's host-side accumulation from the coupler. Port only if profiles show this blocks a GPU stream. |
| #133 | `MOM_wave_interface.F90` — `get_stokessl_lifoxkemper`, `ust_2_u10_coare3p5`, `get_langmuir_number` | 2.44 + 2.32 + 0.84 | Simple scalar math per (i,j). `declare target` + `do concurrent`. |
| #86 | `MOM_spatial_means.F90` — `global_volume_mean`, `global_area_integral`, `global_area_mean`, `array_global_min_max` | 4.20 + 0.26 + 0.17 + 0.03 | Reductions. Use `!$omp target teams distribute ... reduction()`. MPI allreduce still CPU-side — sync the scalar, not the array. |
| #140 | `MOM_tracer_registry.F90` — `post_tracer_transport_diagnostics`, `post_tracer_diagnostics_at_sync`, `postale_tracer_diagnostics` | 1.24 + 0.09 + 0.08 | Diagnostics — only port if flagged on; otherwise accept host-side cost. |
| #127 | `MOM_diagnostics.F90` — `post_transport_diagnostics`, `post_surface_thermo_diags`, `calculate_energy_diagnostics`, `calculate_diagnostic_fields` | 7.30 + 6.03 + 4.99 + 3.38 | Diagnostic-only; sync costs live here. Accept CPU for now. |
| #147 | `MOM_tracer_flow_control.F90` — `call_tracer_set_forcing`, `call_tracer_column_fns`, `call_tracer_surface_state` | 0.11 + 0.06 + 0.07 | Issue body speculates these may need to stay CPU (coupler dispatch) — confirm before spending effort. |

---

## Dependency graph (condensed)

```
Tier 0 (EOS static dispatch, intrinsic math, temperature convert, TFreeze)
 │
 ├── unblocks ──► Tier 1: isopycnal_slopes, neutral_diffusion, ALE/remap/regrid
 │
 ├── unblocks ──► Tier 2: set_viscosity (ML), set_diffusivity, full_convection,
 │ kappa_shear, diagnose_MLD, EPBL (if EOS-dep confirmed)
 │
 └── unblocks ──► Tier 3: thickness_diffuse_full, wave_speed/tdma6,
 density_integrals, PressureForce_Montgomery,
 PressureForce_FV (finishing)

Tier 1 ALE cluster (inter-file):
 ALE ──► remapping ──► regridding
 │ │
 │ └── isopycnal_slopes ◄── neutral_diffusion
 │
 └── ale_remap_tracers ◄── tracer_registry (registry flattening, see §bugs #9)

Tier 2 diabatic chain (sequential in time-step):
 set_viscosity ──► set_diffusivity ──► EPBL / kappa_shear
 └► opacity ──► diabatic_aux/applyboundaryfluxesinout
 └► tracer_diabatic/vertdiff
```

---

## Recommended execution order (concrete)

1. **Fix the 3 remaining confirmed bugs** from §"Suspected bugs" above (items 1–3: CoriolisAdv `id_PV` typo, CoriolisAdv duplicate `rv_x_u`, Montgomery `PFu_bc/PFv_bc` unmapped). One PR each, cheap, restores confidence. Items 4–5 were retracted on verification.
2. **Close the `exit data` gaps** (§bug #6–7). One PR, mechanical.
3. **Tier 0: EOS static dispatch refactor** — coordinate with upstream PR #156 (Roquet). Without this, ~20 Tier 1–3 issues stay blocked or require throwaway workarounds. Expected net wall-clock saving: large, given Roquet density dominates ACCESS-OM3 time.
4. **Tier 0 scalar primitives** (#122, #124, #128) — trivial, can parallelize with step 3 review.
5. **Tier 1 ALE/remap cluster** — user already has context. Port bottom-up: `intrinsic_functions` (#124, done in 4) → `remapping` (#119, user assigned `interpolate_column`) → `regridding` (#97) → `ALE` (#120, user self-assigned). Pair with `isopycnal_slopes` (#75) since neutral diffusion will need it.
6. **Tier 1 neutral_diffusion** (#117). Biggest single EOS-bound hotspot after Roquet itself. Port `neutral_diffusion_calc_coeffs` first, then PPM helpers, then the main `neutral_diffusion` kernel.
7. **Tier 3 heavy hitters by cost** — in order: `mixedlayer_restrat_bodner` (#81, 365s), `density_integrals` (#118, 977s/PLM — Tier 0 must be complete), `wave_speed/tdma6` (#71, 136s).
8. **Tier 2 column physics** — set_viscosity finish (#68, in progress), set_diffusivity (#70), EPBL (#121), kappa_shear (#59, requires refactor per @edoyango comment), then the smaller diabatic_aux / opacity / tidal_mixing items.
9. **Eliminate the wasted roundtrips** (§17, §19 above) — once the CPU routines on the hot path are ported, many of these disappear for free.
10. **Polish**: Tier 4 items, MEKE (#129), diagnostics, wave interface.

---

## Structural items (not a single file — affect the whole port)

These fell outside the "dynamical core + ALE" per-file scope but are worth surfacing because they change the shape of the porting work overall.

### Halo exchange on device (issue #137, `MOM_domain_infra.F90`)

The current pattern in every module is: `!$omp target update from(field)` → `pass_var(field, G%domain)` → `!$omp target update to(field)`. This shows up dozens of times across `MOM_dynamics_split_RK2.F90`, `MOM_barotropic.F90`, `MOM_hor_visc.F90`, etc. It's the single biggest category of unnecessary data movement in the current port — every halo exchange pays 2× full-array transfer on top of the actual MPI cost.

**Options:**
1. **GPU-aware MPI** — call MPI directly on device pointers via `use_device_ptr` clauses. Requires MPI implementation support (OpenMPI+UCX, MPICH+CUDA, Cray MPICH all have it). Eliminates the `update from/to` pairs entirely.
2. **Device-resident halo buffers** — pack/unpack on device into contiguous buffers, only copy buffers to host for MPI. Partial win; still pays the buffer transfer.
3. **Defer** — accept the cost until the per-module port is more complete, then tackle this uniformly. The issue list puts `do_group_pass`, `pass_var_*`, `pass_vector_*` under benchmark config, suggesting the team is tracking this separately.

**Recommendation:** do (1) after Tier 0 lands but before going deep into Tier 2 column physics — the combined reduction in sync points will expose whether what's left is compute-bound or still transfer-bound.

### Diagnostic remap (issue #96, `MOM_diag_remap.F90`)

`do_remap` (120.7s) and `vertically_reintegrate_field` (100.9s) in ACCESS-OM3 are diagnostic-time remapping, not the dynamics-side ALE remap from issue #119. Two separate remap code paths. These share algorithmic DNA with `MOM_remapping.F90` — if the Tier 1 remap port produces reusable `declare target` kernels, they apply here too. 220s combined; worth picking up after #119 lands.

### Surface forcing conversion (issue #125, `mom_surface_forcing_nuopc.F90`)

`convert_iob_to_fluxes` (24.5s) and `convert_iob_to_forces` (4.4s) run once per coupling step. Currently CPU — per §map, the surface forcing fields in `MOM_forcing_type.F90` are intentionally [CPU only] and synced via `update to(U_star)` at line 1233. If the team is OK with the forcing staying CPU-side (it's a coupler boundary), this issue stays deferred. If not, porting makes sense together with #78 (`fluxes_accumulate`, 11.5s).

### Init-time hotspot: `create_depth_list` (issue #116)

`create_depth_list` shows up as **2133.74s** in ACCESS-OM3 — the single largest number in the issue list. It's in `MOM_sum_output.F90`, called once at startup to build a sorted list of cell depths for global energy/mass accounting. **This is not in the per-step loop.** Its cost shows up as wall-clock-at-startup, not throughput. Low-priority for GPU porting (the algorithm is inherently sequential: it's a sort), but worth calling out so the number doesn't mislead planning. If startup time matters, the fix is to make the sort parallel on the host, not to move it to the device.

### Closed issues with load-bearing context

- **#30** (closed): nvfortran 24.9 `-O2` segfaults compiling `MOM_PressureForce_Montgomery.F90`. Check the workaround (flags? refactor?) before re-touching that file — some of the issues in §map (e.g. `PFu_bc` unmapped) may have been hiding behind this compile issue.
- **#31** (closed): nvfortran 25.5 EOS data transfer error. Confirms exactly the polymorphism diagnosis in §EOS: `EOS%type` as abstract `class(...)` re-allocated to concrete types transfers as a zero-size object under 25.5. This is the authoritative evidence that the static-dispatch refactor isn't optional — it's blocking a supported compiler version. Reference it in the Tier 0 PR description.
- **#134, #143, #149**: closed and resolved (safe_alloc, restart, gtracer_flux). No action.

---

## Cross-cutting patterns to reuse

Patterns that recur across the above — worth having a canonical example for each:

- **Column solver port.** Parallelize (j,i), keep k-loop serial inside. Applies to: `ale_plm_edge_values`, `epbl_column`, `kappa_shear_column`, `tdma6`, `make_frazil`, `tridiagts_eulerian`, `tracer_vertdiff`, `streamfn_solver`. Canonical example: find the cleanest one already done and document it.
- **Pure elemental → `declare target`.** Applies to all EOS subclasses, temperature_convert, TFreeze, cuberoot, Langmuir utils.
- **Tridiagonal solve.** Appears in kappa_shear, tracer_vertdiff, diabatic_aux tridiagts, wave_speed tdma6. Serial-k per column; consider Thomas algorithm with careful register pressure.
- **Static-dispatch replacement for `class(...)` polymorphism.** EOS is the main case; `MOM_regridding` grid-kind dispatch is the same shape.
- **Registry/pointer flattening.** Tracer registry (`Reg%Tr(m)%t`) is the canonical blocker; if other modules develop similar patterns (e.g. `OBC%segment(n)%tr_Reg`), apply the same flattening approach.

Closed issue #30 (`-O2` segfault compiling Montgomery pressure force) and #31 (EOS data transfer in Nvidia 25.5) are relevant prior art — check their fixes before re-touching those files.

