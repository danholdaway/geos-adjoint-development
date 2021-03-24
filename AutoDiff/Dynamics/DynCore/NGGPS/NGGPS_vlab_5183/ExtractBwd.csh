#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

cd adm/
cd /gpfsm/dnb31/drholdaw/GEOSadas-5_18_0-vlabfv3pert/src/GEOSagcmPert_GridComp/GEOSdynamicsPert_GridComp/fvdycorepert/model_adm

rm -f AllBWD.txt new.F90 allbwdsubs.F90

#Create list of subroutines that need to be copied later
grep -i "routine" *0 | grep -i "_bwd(" > AllBWD.txt

touch allbwdsubs.F90

foreach line ( "`cat AllBWD.txt`" )

  set file_sen = `echo $line | cut -d ':' -f 1`

  echo $file_sen

  set routine = `echo $line | cut -d ' ' -f 3 | cut -d '(' -f 1`

  echo $routine

  set lsstr = `echo "SUBROUTINE" $routine"("`
  set lestr = `echo "END SUBROUTINE" $routine`

  set ls_sen = `grep -n "$lsstr" $file_sen | cut -d ':' -f 1`
  set le_sen = `grep -n "$lestr" $file_sen | cut -d ':' -f 1`

  awk "NR>=${ls_sen} && NR<=${le_sen}" $file_sen > sub.txt

  cat sub.txt allbwdsubs.F90 > new.F90

  mv new.F90 allbwdsubs.F90

  rm -rf sub.txt


end

rm -f AllBWD.txt


#Done
