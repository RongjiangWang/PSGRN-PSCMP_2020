FORTRAN code for calculating co- and post-seismic deformation in multi-layered viscoelastic half-space based on the viscoelastic-gravitational dislocation theory.

Highlights:

(1) orthonormal propagator algorithm for numerical stability

(2) finite fault model

(3) gravity effect on deformation

(4) output of complete geophysical observables (displacement, stress/strain, tilt, plate rotation, gravity and geoid changes)

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../PSGRN2020/SourceCode

~>gfortran -o psgrn2020 *.f -O3

to get the excutable code pagrn2020.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Okada, Y., Internal deformation due to shear and tensile faults in a half-space, Bull. Seis. Soc. Am., 82, 1018-1040, 1992.

Wang, R., A simple orthonormalization method for the stable and efficient computation of Green's functions, Bull. Seism. Soc. Am., 89, 733-741, 1999.

Wang, R., F. Lorenzo-Martin and F. Roth (2006), PSGRN/PSCMP - a new code for calculating co- and post-seismic deformation, geoid and gravity changes based on the viscoelastic-gravitational dislocation theory, Computers and Geosciences, 32, 527-541. doi:10.1016/j.cageo.2005.08.006
