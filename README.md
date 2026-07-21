# Personalized Elk 9.5.14

This repository contains a research-oriented fork of the [Elk FP-LAPW code](https://elk.sourceforge.io/), based on Elk version 9.5.14.

The fork adds routines for analysing spin textures, orbital- and band-resolved density matrices, modified Fermi surfaces, and momentum- and band-resolved magnetisation densities. The complete Elk source tree is located in `elk-9.5.14/`.

## Research extensions

The custom functionality is available through the following Elk tasks:

| Task | Main routine | Purpose | Main output |
| --- | --- | --- | --- |
| `24` | `spintexture1d` | Calculates the spin texture along a reciprocal-space path using second-variational eigenvectors. | `BAND.OUT`, `BANDLINES.OUT`, `C_EVECSV_INFO.OUT` |
| `26` | `rhonkLOBasis` | Calculates momentum- and band-resolved density matrices in the local-orbital basis. Bands are selected when they remain within ±0.5 Hartree of the Fermi energy along the path. | `BAND.OUT`, `BANDLINES.OUT`, `C_rhonk_INFO.OUT` |
| `103` | `fermisurf` | Writes band-resolved Fermi-surface data on a three-dimensional k-point grid. This version was adapted for research on `Sr2RhO4` and `Sr2RuO4`. | `FERMISURF.OUT`, or spin-resolved Fermi-surface files |
| `901` | `compute_plane_magn_nk_resolved` | Calculates the magnetisation density for individual bands and k-points in a real-space plane. | `MAG2D_ist_<state>_ik_<k-point>.OUT`, `BANDLINES.OUT` |

The fork also supports the `trsym` input variable. Setting it to `.true.` forces the calculated state to be time-reversal symmetric and therefore non-magnetic. This can help when unwanted magnetisation is produced by initial conditions.

The custom implementation is integrated into the normal Elk build. Important additional or modified source files include:

- `elk-9.5.14/src/spintexture1d.f90`
- `elk-9.5.14/src/rhonkLOBasis.f90`
- `elk-9.5.14/src/compute_plane_magn_nk_resolved.f90`
- `elk-9.5.14/src/compute_rhomag_nk_Sebbe.f90`
- `elk-9.5.14/src/rhomagsh_Sebbe.f90`
- `elk-9.5.14/src/vecplot_nk_Sebbe.f90`
- `elk-9.5.14/src/modsebbe.f90`
- modified task dispatch and Fermi-surface code in `elk-9.5.14/src/elk.f90` and `elk-9.5.14/src/fermisurf.f90`

## Requirements

The default configuration in `elk-9.5.14/make.inc` uses:

- a Fortran compiler with MPI support, normally `mpif90`;
- BLAS and LAPACK, configured here with OpenBLAS;
- FFTW3 single- and double-precision libraries;
- OpenMP support.

The default linker flags are:

```text
-lopenblas -llapack -lfftw3 -lfftw3f
```

Edit `elk-9.5.14/make.inc` if your compiler or installed libraries use different names or paths. The file also contains alternative configurations for Intel compilers, MKL, BLIS, serial builds, Libxc, and Wannier90.

## Building

From the repository root:

```bash
cd elk-9.5.14
make all
```

The main executable is created at:

```text
elk-9.5.14/src/elk
```

The build also creates the auxiliary `eos` and `spacegroup` programs in their respective `src` subdirectories.

To remove compiled objects and executables:

```bash
make clean
```

## Running a calculation

Elk reads an input file named `elk.in` from the calculation directory. The required ground-state files must already be present before running post-processing tasks such as tasks `24`, `26`, `103`, or `901`.

For example, an input file can select the spin-texture task with:

```text
tasks
  24

```

The blank line after the task list is required by Elk. Run the calculation from the directory containing `elk.in`:

```bash
/path/to/this/repository/elk-9.5.14/src/elk
```

For an MPI build, launch Elk with the desired number of processes, for example:

```bash
mpirun -np 4 /path/to/this/repository/elk-9.5.14/src/elk
```

To force a non-magnetic, time-reversal-symmetric state, add the following to `elk.in`:

```text
trsym
  .true.
```

The standard Elk input syntax and task descriptions are documented in the source tree and in `elk-9.5.14/README`.

## Testing

The original Elk test targets are available from the `elk-9.5.14` directory:

```bash
cd elk-9.5.14
make test
```

For an MPI build, use:

```bash
make test-mpi
```

Libxc-specific tests are available when Libxc has been enabled in `make.inc`:

```bash
make test-libxc
make test-libxc-mpi
```

## Repository layout

```text
.
├── README.md                 Project-specific documentation
└── elk-9.5.14/
    ├── src/                  Elk and research-modified Fortran sources
    ├── examples/             Example Elk calculations
    ├── tests/                Standard Elk tests
    ├── species/              Atomic species files
    ├── make.inc              Compiler and library configuration
    ├── Makefile              Top-level build and test targets
    ├── README                Original Elk documentation
    └── COPYING               GNU General Public License
```

## License

The bundled Elk code is distributed under the GNU General Public License. See [`elk-9.5.14/COPYING`](elk-9.5.14/COPYING) for the license text.

Research results produced with the custom routines should be checked carefully against appropriate convergence tests and, where possible, against the corresponding standard Elk tasks.
