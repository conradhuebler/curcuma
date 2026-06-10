! Compute tblite's GFN2-D4 ATM three-body energy (and the q=0 two-body baseline).
!
! Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
! Claude Generated 2026. GPL-3.0.
!
! tblite's GFN2 D4 (s9=5.0) splits its energy: the SCC two-body part lives in
! eelec (get_energy at converged charges), while the NON-self-consistent part —
! the q=0 two-body baseline AND the ATM three-body term — is the pre-SCF
! `dispersion` container (get_engrad -> get_dispersion_nonsc). dftd4 evaluates the
! ATM (get_dispersion3) at q=0 C6 (disp.f90:110-116), so the ATM is
! charge-independent.
!
! This tool builds the exact GFN2 model + damping (s6=1,s8=2.7,a1=0.52,a2=5,s9=5,
! cutoffs disp2=50/disp3=25 as in tblite new_d4_dispersion) and reports:
!   E_2body(q=0)  = get_dispersion2 at q=0 C6
!   E_ATM(3body)  = get_dispersion3 at q=0 C6   <-- curcuma's D4Evaluator omits this
! E_2body(q=0)+E_ATM should equal tblite's pre-SCF `dispersion` container
! (ecomp_dispersion in the dump). Curcuma (no ATM) is missing E_ATM.
!
!   dump_dftd4_atm <in.xyz>
program dump_dftd4_atm
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_read, only : read_structure
   use dftd4, only : d4_model, new_d4_model, get_coordination_number
   use dftd4_model, only : d4_ref
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_damping_rational, only : rational_damping_param
   implicit none

   type(structure_type) :: mol
   type(d4_model) :: model
   type(rational_damping_param) :: param
   type(error_type), allocatable :: error
   type(realspace_cutoff) :: cutoff
   character(len=4096) :: arg
   integer :: nat, mref
   real(wp), allocatable :: lattr(:, :), cn(:), q(:), gwvec(:, :), c6(:, :)
   real(wp), allocatable :: e2(:), e3(:)
   real(wp) :: e2body, eatm

   if (command_argument_count() < 1) then
      write(0, '(a)') "usage: dump_dftd4_atm <in.xyz>"
      stop 2
   end if
   call get_command_argument(1, arg)
   call read_structure(mol, trim(arg), error)
   if (allocated(error)) then
      write(0, '(a)') "read_structure failed: "//error%message
      stop 3
   end if
   nat = mol%nat

   call new_d4_model(model, mol, ref=d4_ref%gfn2)
   ! GFN2 damping (tblite xtb/gfn2.f90:54) + cutoffs (disp/d4.f90:82).
   param = rational_damping_param(s6=1.0_wp, s8=2.7_wp, s9=5.0_wp, a1=0.52_wp, a2=5.0_wp)
   cutoff = realspace_cutoff()
   cutoff%disp2 = 50.0_wp
   cutoff%disp3 = 25.0_wp
   mref = maxval(model%ref)

   ! CN (tblite uses model%rcov/en).
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   allocate(cn(nat))
   call get_coordination_number(mol, lattr, cutoff%cn, model%rcov, model%en, cn)

   ! q=0 reference C6 (dftd4 computes BOTH the ATM and — in this tool — the
   ! baseline two-body at q=0; the SCC two-body correction is separate).
   allocate(q(nat)); q(:) = 0.0_wp
   allocate(gwvec(mref, nat))
   call model%weight_references(mol, cn, q, gwvec)
   allocate(c6(nat, nat))
   call model%get_atomic_c6(mol, gwvec, c6=c6)

   ! Two-body at q=0.
   allocate(e2(nat)); e2(:) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, model%r4r2, c6, energy=e2)
   e2body = sum(e2)

   ! ATM three-body.
   allocate(e3(nat)); e3(:) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, model%r4r2, c6, energy=e3)
   eatm = sum(e3)

   write(*, '(a,i0)')        "nat            = ", nat
   write(*, '(a,es20.12,a)') "E_2body(q=0)   = ", e2body, " Eh"
   write(*, '(a,es20.12,a)') "E_ATM(3body)   = ", eatm,   " Eh   <-- curcuma omits this"
   write(*, '(a,es20.12,a)') "E_2body+E_ATM  = ", e2body+eatm, " Eh (= tblite pre-SCF dispersion container)"

end program dump_dftd4_atm
