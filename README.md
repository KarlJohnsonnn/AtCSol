# AtCSol
We present version 1.0 of the Atmospheric Chemistry Solver - AtCSol -  a fully modularised modern Fortran software package for the computational simulation of chemical kinetic systems derived from atmospheric and combustion chemistry. The purpose of AtCSol is the investigation of various mathematical methods applied to chemistry problems. Version 1.0 of AtCSol contains of a box-model framework, where a mono-disperse distribution of cloud droplets is assumed, several stiff numerical integrators with different measures of local errors, an analytical Jacobian matrix approach and two different ways to solve the upcoming linear equation systems. Next to atmospheric multiphase chemistry problems AtCSol is capable of analysing detailed combustion mechanisms as well, where the time span of interest is several orders of magnitude smaller then atmospheric systems. This is the first public available version of AtCSol and is published under the GNU General Public Licence.


Installation:

  1.  run the shell script in the SOLVER folder:  ./SOLVER/compile_odeSolver.sh
  
  2.  open M_DEF_Mac and set the corret paths to NetCDF and LAPack
  
  3.  run: make
  
      3.1 make AtCSol     (for optimized version)
  
      3.2 make AtCSol_dbg (for debugging version)
      
  4.  run: make test      (for several test scenarios)
  
  5. excecute a specific mechanism:   ./AtCSol.exe RUN/[mechanism].run
