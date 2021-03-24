#!/bin/csh

reset


$TAPENADE_HOME/bin/tapenade -tangent -head radcouple -outvars "RAD_CF RAD_QL RAD_QI RAD_RL RAD_RI" -vars "TE CF AF QCLLS QCILS QCLAN QCIAN" radcouple.f90 -output radcouple -html

$TAPENADE_HOME/bin/tapenade -reverse -head radcouple -outvars "RAD_CF RAD_QL RAD_QI RAD_RL RAD_RI" -vars "TE CF AF QCLLS QCILS QCLAN QCIAN" radcouple.f90 -output radcouple -html


cat combine/top radcouple.f90 radcouple_d.f90 radcouple_b.f90 combine/foot > cloudradcouple.F90

rm -f radcouple_d.F90 radcouple_b.F90

cp cloudradcouple.F90 $radpertdir

