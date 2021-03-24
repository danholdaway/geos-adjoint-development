#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -forward  -inputlanguage fortran95 -outputlanguage fortran95 -head tracer_2d_1L_dh -outvars "q dp0 q3" -vars "q dp0 mfx mfy cx cy q3" tracer_2d.f90 -output tracer_2d_1L_dh -html

$TAPENADE_HOME/bin/tapenade -backward -inputlanguage fortran95 -outputlanguage fortran95 -head tracer_2d_1L_dh -outvars "q dp0 q3" -vars "q dp0 mfx mfy cx cy q3" tracer_2d.f90 -output tracer_2d_1L_dh -html

# Remove the html logs
# --------------------
rm -r tapenadehtml

# Remove the messages
# -------------------
#rm cloud_d.msg cloud_b.msg

# Remove tmp files
# ----------------
rm *~


