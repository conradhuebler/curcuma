! Dump the GFN2-D4 reference C6 matrix from dftd4's new_d4_model.
!
! Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
! Claude Generated 2026. GPL-3.0.
!
! Reads an xyz file via mctc-io, constructs the *exact* d4_model that
! tblite's GFN2 dispersion container builds (new_d4_model(..., ref=d4_ref%gfn2),
! disp/d4.f90:80), and writes the reference C6 matrix model%c6(iref,jref,isp,jsp)
! to JSON. This is the geometry/charge-INDEPENDENT reference C6 (set once during
! model construction from the zeta-scaled aiw via Casimir-Polder), the suspected
! source of curcuma's ~5% C-path residual (see docs/GFN2_D4_STATUS.md).
!
! The companion C++ tool diag_curcuma_d4_c6 reads this JSON and diffs against
! curcuma's D4ParameterGenerator::getC6FlatCache() per (Z_i, Z_j, iref, jref).
!
!   dump_dftd4_c6 <in.xyz> [out.json]   (stdout if out.json omitted)
!
! JSON layout:
!   { "nid": <n species>, "mref": <max refs>,
!     "species_z": [Z per species],
!     "nref": [n refs per species],
!     "c6": [ flat nid*nid*mref*mref, C-order (((isp*nid+jsp)*mref+iref)*mref+jref) ] }
program dump_dftd4_c6
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_read, only : read_structure
   use dftd4, only : d4_model, new_d4_model
   use dftd4_model, only : d4_ref
   implicit none

   type(structure_type) :: mol
   type(d4_model) :: model
   type(error_type), allocatable :: error
   character(len=:), allocatable :: xyz_path, out_path
   character(len=4096) :: arg
   integer :: nargs, unit, isp, jsp, iref, jref, mref, nid
   logical :: first

   nargs = command_argument_count()
   if (nargs < 1) then
      write(0, '(a)') "usage: dump_dftd4_c6 <in.xyz> [out.json]"
      stop 2
   end if
   call get_command_argument(1, arg)
   xyz_path = trim(arg)
   if (nargs >= 2) then
      call get_command_argument(2, arg)
      out_path = trim(arg)
   end if

   call read_structure(mol, xyz_path, error)
   if (allocated(error)) then
      write(0, '(a)') "read_structure failed: "//error%message
      stop 3
   end if

   ! Identical to tblite/src/tblite/disp/d4.f90:80 (new_d4_dispersion).
   call new_d4_model(model, mol, ref=d4_ref%gfn2)

   nid = mol%nid
   mref = maxval(model%ref)

   if (allocated(out_path)) then
      ! Large record length: the c6 array can be thousands of chars on one line.
      open(newunit=unit, file=out_path, action='write', status='replace', recl=1048576)
   else
      unit = 6
   end if

   write(unit, '(a)')      "{"
   write(unit, '(a)')      '  "source": "dftd4 new_d4_model ref=gfn2 model%c6",'
   write(unit, '(a,i0,a)') '  "nid": ', nid, ','
   write(unit, '(a,i0,a)') '  "mref": ', mref, ','

   ! species Z list (1-based species index -> Z)
   write(unit, '(a)', advance='no') '  "species_z": ['
   do isp = 1, nid
      if (isp > 1) write(unit, '(a)', advance='no') ', '
      write(unit, '(i0)', advance='no') mol%num(isp)
   end do
   write(unit, '(a)') '],'

   ! number of reference systems per species
   write(unit, '(a)', advance='no') '  "nref": ['
   do isp = 1, nid
      if (isp > 1) write(unit, '(a)', advance='no') ', '
      write(unit, '(i0)', advance='no') model%ref(isp)
   end do
   write(unit, '(a)') '],'

   ! reference C6 matrix, C-order flat:
   !   idx = (((isp-1)*nid + (jsp-1))*mref + (iref-1))*mref + (jref-1)
   ! entries beyond nref(isp)/nref(jsp) are zero (model%c6 allocated mref x mref).
   ! One value per line to avoid any formatted-record-length limit on the
   ! (potentially thousands of entries) flat c6 array. Entries with iref/jref
   ! beyond an element's reference count are written as 0 — model%c6 is allocated
   ! mref x mref but only filled up to ref(isp), so those slots are uninitialized
   ! memory (denormals). The differ never reads them, but they must be valid JSON.
   ! es24.16e3 forces a 3-digit exponent WITH the 'E' (plain es drops 'E' for
   ! |exp|>=100, producing invalid JSON like 6.0...-154).
   write(unit, '(a)') '  "c6": ['
   first = .true.
   do isp = 1, nid
      do jsp = 1, nid
         do iref = 1, mref
            do jref = 1, mref
               if (.not. first) write(unit, '(a)') ','
               first = .false.
               if (iref <= model%ref(isp) .and. jref <= model%ref(jsp)) then
                  write(unit, '(a,es24.16e3)', advance='no') '    ', model%c6(iref, jref, isp, jsp)
               else
                  write(unit, '(a)', advance='no') '    0.0'
               end if
            end do
         end do
      end do
   end do
   write(unit, '(a)') ''
   write(unit, '(a)') '  ]'
   write(unit, '(a)') "}"

   if (allocated(out_path)) close(unit)

end program dump_dftd4_c6
