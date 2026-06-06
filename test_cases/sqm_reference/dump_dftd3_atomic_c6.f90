! Dump simple-dftd3's weighted per-atom C6 matrix (get_atomic_c6) for a molecule.
!
! Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
! Claude Generated 2026. GPL-3.0.
!
! Reproduces exactly what tblite's GFN1 D3 dispersion does (s-dftd3 disp.f90:171-178):
!   new_d3_model(model, mol)
!   get_coordination_number(mol, lattr, cutoff%cn, model%rcov, cn)
!   model%weight_references(mol, cn, gwvec)
!   model%get_atomic_c6(mol, gwvec, c6=c6)
! Writes the per-atom CN it used and c6(nat,nat) to JSON, so the C++ differ can
! compare curcuma's interpolated D3 C6 (D3ParameterGenerator) against s-dftd3 --
! the term WP2 isolated as the GFN1 dispersion residual (every other D3 term --
! CN, C8/C6, BJ R0, energy formula -- already verified bit-matching s-dftd3).
!
!   dump_dftd3_atomic_c6 <in.xyz> [out.json]
!
! JSON: { "nat":N, "cn":[...], "c6":[ N*N row-major c6(i,j) ] }
program dump_dftd3_atomic_c6
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_read, only : read_structure
   use dftd3, only : d3_model, new_d3_model, get_coordination_number, &
      & realspace_cutoff, get_lattice_points
   implicit none

   type(structure_type) :: mol
   type(d3_model) :: model
   type(error_type), allocatable :: error
   type(realspace_cutoff) :: cutoff
   character(len=:), allocatable :: xyz_path, out_path
   character(len=4096) :: arg
   integer :: nargs, unit, i, j, mref, nat
   real(wp), allocatable :: lattr(:, :), cn(:), gwvec(:, :), c6(:, :)
   logical :: first

   nargs = command_argument_count()
   if (nargs < 1) then
      write(0, '(a)') "usage: dump_dftd3_atomic_c6 <in.xyz> [out.json]"
      stop 2
   end if
   call get_command_argument(1, arg); xyz_path = trim(arg)
   if (nargs >= 2) then
      call get_command_argument(2, arg); out_path = trim(arg)
   end if

   call read_structure(mol, xyz_path, error)
   if (allocated(error)) then
      write(0, '(a)') "read_structure failed: "//error%message
      stop 3
   end if
   nat = mol%nat

   call new_d3_model(model, mol)
   mref = maxval(model%ref)

   ! CN exactly as s-dftd3's energy path (disp.f90:171-172).
   cutoff = realspace_cutoff()
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   allocate(cn(nat))
   call get_coordination_number(mol, lattr, cutoff%cn, model%rcov, cn)

   allocate(gwvec(mref, nat))
   call model%weight_references(mol, cn, gwvec)

   allocate(c6(nat, nat))
   call model%get_atomic_c6(mol, gwvec, c6=c6)

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
   write(unit, '(a)') '  ],'

   ! Per-species reference blocks: Z, nref, reference CN, and the self C6 block
   ! model%c6(ri,rj,isp,isp). Ground truth for diffing curcuma's loaded D3 tables.
   write(unit, '(a)') '  "refdata": {'
   block
      integer :: isp, ri, rj, nref
      logical :: firstsp
      firstsp = .true.
      do isp = 1, mol%nid
         if (.not. firstsp) write(unit, '(a)') ','
         firstsp = .false.
         nref = model%ref(isp)
         write(unit, '(a,i0,a)') '    "', mol%num(isp), '": {'
         write(unit, '(a,i0,a)') '      "nref": ', nref, ','
         write(unit, '(a)', advance='no') '      "refcn": ['
         do ri = 1, nref
            if (ri > 1) write(unit, '(a)', advance='no') ', '
            write(unit, '(es20.12e3)', advance='no') model%cn(ri, isp)
         end do
         write(unit, '(a)') '],'
         write(unit, '(a)') '      "c6self": ['
         do ri = 1, nref
            write(unit, '(a)', advance='no') '        ['
            do rj = 1, nref
               if (rj > 1) write(unit, '(a)', advance='no') ', '
               write(unit, '(es22.14e3)', advance='no') model%c6(ri, rj, isp, isp)
            end do
            if (ri < nref) then
               write(unit, '(a)') '],'
            else
               write(unit, '(a)') ']'
            end if
         end do
         write(unit, '(a)', advance='no') '      ]'
         write(unit, '(a)', advance='no') '    }'
      end do
      write(unit, '(a)') ''
   end block
   write(unit, '(a)') '  }'
   write(unit, '(a)') "}"

   if (allocated(out_path)) close(unit)

end program dump_dftd3_atomic_c6
