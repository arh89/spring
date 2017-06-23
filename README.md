# Spring

Spring is a software package designed to be a simple testbed for atomistic
simulation methods - with particular emphasis on path integral molecular
dynamics based methods.

The original focus of this software was to provide an implementation of the
partially adiabatic centroid molecular dynamics method - where code simplicity
and reusability (rather than efficiency) should take precedence in any design
decisions.

## Known issues

The package contains a very rough (half finished) implementation of a lattice
dynamics code, used to calculate the phonons of a material (through the method
of finite displacements). This is known to contain bugs and should not be used.

The potentials contained in the `potentials` directory are suspected to be
incorrect and/or not integrated with the rest of the code.

Due to the current naming convention of modules, this code will not compile
with Intel Fortran compilers (only tested and known to work with gfortran).

## Requirements

gfortran, LAPACK (OpenBLAS preferred), FFTW3

## Usage

No documentation is currently provided. Please see `inputs/input_keywords` for a
list of valid keywords to the input file.

Please see `inputs/dft_pes.in` for an example input file.
