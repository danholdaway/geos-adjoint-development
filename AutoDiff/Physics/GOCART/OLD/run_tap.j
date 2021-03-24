#!/bin/csh

reset

#$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "emissions" -vars "u10m v10m" DustEmissionGOCART.f90 -output DustEmissionGOCART_Pert -tgtvarname _tl -tgtfuncname _tlm

#$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "emissions" -vars "u10m v10m" DustEmissionGOCART.f90 -output DustEmissionGOCART_Pert -adjvarname _ad -adjfuncname _adm


#$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "AERO" -vars "tmpu RH rhoa AERO" Chem_Settling.f90 -output Chem_Settling_Pert -tgtvarname _tl -tgtfuncname _tlm

#$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "AERO" -vars "tmpu RH rhoa AERO" Chem_Settling.f90 -output Chem_Settling_Pert -adjvarname _ad -adjfuncname _adm


$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "drydepf" -vars "tmpu u10m v10m rhoa" DryDepositionGOCART.f90 -output DryDepositionGOCART_Pert -tgtvarname _tl -tgtfuncname _tlm

$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "drydepf" -vars "tmpu u10m v10m rhoa" DryDepositionGOCART.f90 -output DryDepositionGOCART_Pert -adjvarname _ad -adjfuncname _adm

#---------------------------------------

#$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "rho" -vars "PT QV" airdens.f90 -output airdens -tgtvarname _tl -tgtfuncname _tlm

#$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran90 -outputlanguage fortran90 -head dummy -outvars "rho" -vars "PT QV" airdens.f90 -output airdens -adjvarname P -adjfuncname _adm

