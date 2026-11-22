# CPCM-X solvent files for cli_sqm_21

`water.cosmo`, `hexane.cosmo` and `acetonitrile.cosmo` are copied unchanged from `DB/xtb/` of CPCM-X 1.1.0
(https://github.com/grimme-lab/CPCM-X, commit e7f894c, M. Stahn, S. Ehlert, S. Grimme), licensed under the GNU Lesser
General Public License v3.0 or later. They are the conductor-COSMO sigma segments of the solvents computed with GFN2-xTB.
Only these three files are kept here (188 KB); the full database (88 solvents, 8.9 MB of .cosmo files) is not part of
curcuma: point `-xtb.cpcmx_db` or `CURCUMA_CPCMX_DB` to a checkout of `CPCM-X/DB/xtb` (see docs/SOLVATION.md).
