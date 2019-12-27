##TAM 570 Computational Fluid Mechanics (Fall 2019)

A repository containing code developed for TAM 570 at UIUC, with my partners Anand Radhakrishnan and Dongchen Wang. Several utility codes provided by Prof. Paul Fischer. The course used spectral methods to solve fluid dynamics-related problems, including a Laminar Navier-Stokes solver.

All files were developed on and worked with MATLAB 2019.

**Important:** With MATLAB 2019, SPD matrices must be explicitly enforced i.e., set `A = 0.5*(A+A')` if `A` is SPD, even if the routine that generates A is supposed to generate an SPD matrix. Can cause a lot of heartaches if this isn't done!
