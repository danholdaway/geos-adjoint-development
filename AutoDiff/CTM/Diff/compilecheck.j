#!/bin/csh
 
reset

source /gpfsm/dnb31/drholdaw/Heracles-5_4_p3_CTM-r1/Linux/bin/g5_modules

#Clean up
rm -f *.mod *.o

ifort -check -c ESMF.F90

ifort -check -c GmiArrayBundlePointer_mod.F90

ifort -check -c updateDiffusion_mod.F90


#Clean up
rm -f *.mod *.o
