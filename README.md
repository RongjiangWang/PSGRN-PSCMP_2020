FORTRAN code for calculating co- and post-seismic deformation in multi-layered viscoelastic half-space based on the viscoelastic-gravitational dislocation theory.

Highlights:

(1) orthonormal propagator algorithm for numerical stability

(2) finite fault model

(3) azimuthal equidistant projection between spherical and local cartesian coordinates

(4) gravity effect on deformation

(5) output of complete geophysical observables (displacement, stress/strain, tilt, plate rotation, gravity and geoid changes)

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../PSGRN2020/SourceCode

~>gfortran -o psgrn2020 *.f -O3

to get the excutable code psgrn2020.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Okada, Y., Internal deformation due to shear and tensile faults in a half-space, Bull. Seis. Soc. Am., 82, 1018-1040, 1992.

Wang, R., A simple orthonormalization method for the stable and efficient computation of Green's functions, Bull. Seism. Soc. Am., 89, 733-741, 1999.

Wang, R., F. Lorenzo-Martin and F. Roth (2006), PSGRN/PSCMP - a new code for calculating co- and post-seismic deformation, geoid and gravity changes based on the viscoelastic-gravitational dislocation theory, Computers and Geosciences, 32, 527-541. doi:10.1016/j.cageo.2005.08.006.

Wang, R. (2005), On the singularity problem of the elastic-gravitational dislocation theory applied to plane-Earth model, Geophysical Research Letters, 32, L06307, doi:10.1029/2003GL019358.

Wang, R. (2005), The dislocation theory: a consistent way for including the gravity effect in (vis-co)elastic plane-earth models, Geophysical Journal International, 161, 191-196.

Wang, R. (2007), Erratum to “The dislocation theory: a consistent way for including the gravity effect in (visco)elastic plane-earth models”, Geophysical Journal International, 170, 857.

---------------------------------------------------------
Last update on Nov 24, 2025, in Zhuhai, by Rongjiang Wang
