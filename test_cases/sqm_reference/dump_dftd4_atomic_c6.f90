! Dump tblite's WEIGHTED per-atom C6 matrix (get_atomic_c6) at given charges.
!
! Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
! Claude Generated 2026. GPL-3.0.
!
! Reproduces exactly what tblite's GFN2 dispersion get_engrad does
! (disp/d4.f90:219-224):
!   new_d4_model(model, mol, ref=d4_ref%gfn2)
!   get_coordination_number(mol, lattr, cutoff%cn, model%rcov, model%en, cn)
!   model%weight_references(mol, cn, qat, gwvec, ...)
!   model%get_atomic_c6(mol, gwvec, ..., c6)
! where qat are tblite's converged Mulliken charges (read from <charges.txt>,
! nat whitespace-separated doubles in input order). Writes c6(nat,nat) + the
! per-atom CN it used to JSON, so the C++ differ can compare curcuma's
! weightedC6Gfn2 AND curcuma's D4 CN against tblite's.
!
!   dump_dftd4_atomic_c6 <in.xyz> <charges.txt> [out.json]
!
! JSON: { "nat":N, "cn":[...], "qat":[...], "c6":[ N*N row-major c6(i,j) ] }
program dump_dftd4_atomic_c6
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_read, only : read_structure
   use dftd4, only : d4_model, new_d4_model, get_coordination_number
   use dftd4_model, only : d4_ref
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   implicit none

   type(structure_type) :: mol
   type(d4_model) :: model
   type(error_type), allocatable :: error
   type(realspace_cutoff) :: cutoff
   character(len=:), allocatable :: xyz_path, q_path, out_path
   character(len=4096) :: arg
   integer :: nargs, unit, qunit, i, j, mref, nat, ios
   real(wp), allocatable :: lattr(:, :), cn(:), gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :), qat(:)
   logical :: first

   nargs = command_argument_count()
   if (nargs < 2) then
      write(0, '(a)') "usage: dump_dftd4_atomic_c6 <in.xyz> <charges.txt> [out.json]"
      stop 2
   end if
   call get_command_argument(1, arg); xyz_path = trim(arg)
   call get_command_argument(2, arg); q_path = trim(arg)
   if (nargs >= 3) then
      call get_command_argument(3, arg); out_path = trim(arg)
   end if

   call read_structure(mol, xyz_path, error)
   if (allocated(error)) then
      write(0, '(a)') "read_structure failed: "//error%message
      stop 3
   end if
   nat = mol%nat

   ! Read converged charges (nat doubles, input/atom order).
   allocate(qat(nat))
   open(newunit=qunit, file=q_path, action='read', status='old')
   do i = 1, nat
      read(qunit, *, iostat=ios) qat(i)
      if (ios /= 0) then
         write(0, '(a,i0)') "charges.txt read error at atom ", i
         stop 4
      end if
   end do
   close(qunit)

   call new_d4_model(model, mol, ref=d4_ref%gfn2)
   mref = maxval(model%ref)

   ! CN exactly as tblite's dispersion update (disp/d4.f90:105-107).
   cutoff = realspace_cutoff()
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   allocate(cn(nat))
   call get_coordination_number(mol, lattr, cutoff%cn, model%rcov, model%en, cn)

   allocate(gwvec(mref, nat), gwdcn(mref, nat), gwdq(mref, nat))
   call model%weight_references(mol, cn, qat, gwvec, gwdcn, gwdq)

   ! Diagnostic: dump gwvec for the first atom of each element (stderr), to
   ! compare per-reference weights against curcuma's buildRefW.
   block
      integer :: a, r
      do a = 1, min(nat, 40)
         write(0, '(a,i0,a,i0,a)', advance='no') "GWVEC atom ", a, " Z=", mol%num(mol%id(a)), " :"
         do r = 1, model%ref(mol%id(a))
            write(0, '(1x,es16.8)', advance='no') gwvec(r, a)
         end do
         write(0, '(a)') ""
      end do
   end block

   allocate(c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat))
   call model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   if (allocated(out_path)) then
      open(newunit=unit, file=out_path, action='write', status='replace', recl=1048576)
   else
      unit = 6
   end if

   write(unit, '(a)')      "{"
   write(unit, '(a,i0,a)') '  "nat": ', nat, ','

   write(unit, '(a)', advance='no') '  "cn": ['
   do i = 1, nat
      if (i > 1) write(unit, '(a)', advance='no') ', '
      write(unit, '(es22.14e3)', advance='no') cn(i)
   end do
   write(unit, '(a)') '],'

   write(unit, '(a)', advance='no') '  "qat": ['
   do i = 1, nat
      if (i > 1) write(unit, '(a)', advance='no') ', '
      write(unit, '(es22.14e3)', advance='no') qat(i)
   end do
   write(unit, '(a)') '],'

   ! c6(i,j) row-major: idx = (i-1)*nat + (j-1). One value per line.
   write(unit, '(a)') '  "c6": ['
   first = .true.
   do i = 1, nat
      do j = 1, nat
         if (.not. first) write(unit, '(a)') ','
         first = .false.
         write(unit, '(a,es24.16e3)', advance='no') '    ', c6(i, j)
      end do
   end do
   write(unit, '(a)') ''
   write(unit, '(a)') '  ]'
   write(unit, '(a)') "}"

   if (allocated(out_path)) close(unit)

end program dump_dftd4_atomic_c6
