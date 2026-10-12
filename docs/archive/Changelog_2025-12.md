# Changelog

## [0.0.190] - 2025-12-13
### Fixed
- GFNFF torsion generation:
  - Fixed unit mismatch in bond detection (Angstrom vs Bohr units)
  - Fixed duplicate torsion detection algorithm to properly compare all four atoms
  - Corrected torsion count for ethane from 36 (incorrect) to 9 (correct)
  - Fixed neighbor list generation for accurate torsion topology

## [0.0.189] - 2025-12-12
### Added
- Initial GFNFF implementation with bond, angle, and torsion parameters
- Support for JSON parameter files
- Basic energy calculations for small molecules

[Unreleased]: https://github.com/conradhuebler/curcuma/compare/v0.0.190...HEAD
[0.0.190]: https://github.com/conradhuebler/curcuma/compare/v0.0.189...v0.0.190
[0.0.189]: https://github.com/conradhuebler/curcuma/releases/tag/v0.0.189