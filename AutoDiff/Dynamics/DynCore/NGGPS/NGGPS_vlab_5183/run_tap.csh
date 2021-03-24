#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


##Generate tangent linear code
##----------------------------
#@ start_time = `date +%s`
#
#rm -rf tlm/
#mkdir tlm/
#rm -rf tlm/tapenadehtml
#
#$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "main_mod.run (u v pt delp q w delz)/(u v pt delp q w delz)" main.F90 fv_dynamics.F90 a2b_edge.F90 boundary.F90 dyn_core.F90 fv_arrays.F90 fv_arrays_nlm.F90 fv_cmp.F90 fv_control.F90 fv_fill.F90 fv_grid_utils.F90 fv_mapz.F90 fv_nesting.F90 fv_sg.F90 fv_tracer2d.F90 nh_core.F90 nh_utils.F90 sw_core.F90 tp_core.F90 MAPL/*.F90 tools/*.F90 fms/*.F90
#
#cd tlm/
#
##Rename file ending to match the model
#rename mod_diff.f90 tlm.F90 *mod_diff.f90
#
##Remove the files we dont need
#rm -f mpp* init_hydro_tlm.F90 fv_arrays_tlm.F90 fv_arrays_nlm_tlm.F90 fms_io_tlm.F90
#
###Perform uncommenting
#sed -i -- 's/\!UNCOM/      /g' *_tlm.F90
#
##Replace and rename _FB schemes, where subroutine
#grep -i '_FB_TLM(' *0 | grep -i 'SUBROUTINE' > start.txt
#grep -i '_FB_TLM' *0 | grep -i 'END SUBROUTINE' > final.txt
#paste start.txt final.txt > LinesToDelete.txt
#rm -f start.txt final.txt
#foreach line ( "`cat LinesToDelete.txt`" )
# set file = `echo $line | cut -d ':' -f 1`
# set sub_tl_ren = `echo $line | cut -d ' ' -f 3 | cut -d '(' -f 1`
# set sub_fb = `echo $sub_tl_ren | rev | cut -c 5- | rev`
# set sub_nl = `echo $sub_tl_ren | rev | cut -c 8- | rev`
# set sub_tl_del = `echo ${sub_nl}_TLM`
# set sub_nl_lc = `echo $sub_nl | tr 'A-Z' 'a-z'`
# set sub_fb_lc = `echo $sub_fb | tr 'A-Z' 'a-z'`
# #First delete the TLM of the full scheme
# set s1 = `echo Differentiation of ${sub_nl} in`
# set s2 = `echo END SUBROUTINE ${sub_tl_del}`
# set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
# set le = `grep -in "$s2" $file | cut -d ':' -f 1`
# awk "NR<${ls}" $file > a.txt
# awk "NR>${le}" $file > b.txt
# cat a.txt b.txt > ${file}
# rm -f a.txt b.txt
# #Now rename the _FB_TLM scheme to just _TLM
# set r1 = `echo SUBROUTINE ${sub_tl_ren}`
# set r2 = `echo SUBROUTINE ${sub_tl_del}`
# sed -i -- "s/${r1}/${r2}/g" $file
# #Now delete the nonlinear _FB scheme
# set s1 = `echo SUBROUTINE ${sub_fb}\(`
# set s2 = `echo END SUBROUTINE ${sub_fb}`
# set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
# set le = `grep -in "$s2" $file | cut -d ':' -f 1`
# awk "NR<${ls}" $file > a.txt
# awk "NR>${le}" $file > b.txt
# cat a.txt b.txt > ${file}
# rm -f a.txt b.txt
# #Now rename the calls to _FB_TLM scheme to call just _TLM, all files
# set r1 = `echo CALL ${sub_tl_ren}\(`
# set r2 = `echo CALL ${sub_tl_del}\(`
# sed -i -- "s/${r1}/${r2}/g" *0
# #Now rename the calls to _FB nonlinear scheme to call just nl scheme, all files
# set r1 = `echo CALL ${sub_fb}\(`
# set r2 = `echo CALL ${sub_nl}\(`
# sed -i -- "s/${r1}/${r2}/g" *0
# #Rename lower case parts, uses and Differentiation of...
# sed -i -- "s/${sub_fb_lc}/${sub_nl_lc}/g" *0
#end
#rm -f LinesToDelete.txt
#
##Replace and rename _FB schemes, where function
#grep -i '_FB_TLM(' *0 | grep -i 'FUNCTION' > start.txt
#grep -i '_FB_TLM' *0 | grep -i 'END FUNCTION' > final.txt
#paste start.txt final.txt > LinesToDelete.txt
#rm -f start.txt final.txt
#foreach line ( "`cat LinesToDelete.txt`" )
# set file = `echo $line | cut -d ':' -f 1`
# set sub_tl_ren = `echo $line | cut -d ' ' -f 4 | cut -d '(' -f 1`
# set sub_fb = `echo $sub_tl_ren | rev | cut -c 5- | rev`
# set sub_nl = `echo $sub_tl_ren | rev | cut -c 8- | rev`
# set sub_tl_del = `echo ${sub_nl}_TLM`
# set sub_nl_lc = `echo $sub_nl | tr 'A-Z' 'a-z'`
# set sub_fb_lc = `echo $sub_fb | tr 'A-Z' 'a-z'`
# #First delete the TLM of the full scheme
# set s1 = `echo Differentiation of ${sub_nl} in`
# set s2 = `echo END FUNCTION ${sub_tl_del}`
# set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
# set le = `grep -in "$s2" $file | cut -d ':' -f 1`
# awk "NR<${ls}" $file > a.txt
# awk "NR>${le}" $file > b.txt
# cat a.txt b.txt > ${file}
# rm -f a.txt b.txt
# #Now rename the _FB_TLM scheme to just _TLM
# set r1 = `echo FUNCTION ${sub_tl_ren}`
# set r2 = `echo FUNCTION ${sub_tl_del}`
# sed -i -- "s/${r1}/${r2}/g" $file
# #Now delete the nonlinear _FB scheme
# set s1 = `echo FUNCTION ${sub_fb}\(`
# set s2 = `echo END FUNCTION ${sub_fb}`
# set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
# set le = `grep -in "$s2" $file | cut -d ':' -f 1`
# awk "NR<${ls}" $file > a.txt
# awk "NR>${le}" $file > b.txt
# cat a.txt b.txt > ${file}
# rm -f a.txt b.txt
# #Now rename the calls to _FB_TLM scheme to call just _TLM, all files
# set r1 = `echo = ${sub_tl_ren}\(`
# set r2 = `echo = ${sub_tl_del}\(`
# sed -i -- "s/${r1}/${r2}/g" *0
# #Now rename the calls to _FB nonlinear scheme to call just nl scheme, all files
# set r1 = `echo = ${sub_fb}\(`
# set r2 = `echo = ${sub_nl}\(`
# sed -i -- "s/${r1}/${r2}/g" *0
# #Rename lower case parts, uses and Differentiation of...
# sed -i -- "s/${sub_fb_lc}/${sub_nl_lc}/g" *0
#end
#rm -f LinesToDelete.txt
#
#cd ../
#
#@ end_time = `date +%s`
#@ diff = $end_time - $start_time
#echo "TLM generation took $diff seconds"

#Generate adjoint code
#---------------------
@ start_time = `date +%s`

#Generate the adjoint with all the adm/fwd/bwd routines we need
#--------------------------------------------------------------

rm -rf adm/
mkdir adm/

#$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "main_mod.run (u v pt delp q w delz)/(u v pt delp q w delz)" main.F90 fv_dynamics.F90 a2b_edge.F90 boundary.F90 dyn_core.F90 fv_arrays.F90 fv_arrays_nlm.F90 fv_cmp.F90 fv_control.F90 fv_fill.F90 fv_grid_utils.F90 fv_mapz.F90 fv_nesting.F90 fv_sg.F90 fv_tracer2d.F90 nh_core.F90 nh_utils.F90 sw_core.F90 tp_core.F90 MAPL/*.F90 tools/*.F90 fms/*.F90 -nocheckpoint "a2b_edge_mod.a2b_ord4_fb a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar_fb fv_mapz_mod.map1_ppm_fb fv_mapz_mod.mapn_tracer_fb fv_mapz_mod.map1_q2_fb fv_mapz_mod.remap_2d fv_mapz_mod.scalar_profile_fb fv_mapz_mod.cs_profile_fb fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested fv_sg_mod.fv_subgrid_z main_mod.compute_pressures main_mod.run nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u_fb sw_core_mod.ytp_v_fb sw_core_mod.compute_divergence_damping_fb sw_core_mod.smag_corner_fb tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d_fb tp_core_mod.copy_corners_fb tp_core_mod.xppm_fb tp_core_mod.yppm_fb tp_core_mod.deln_flux_fb a2b_edge_mod.extrap_corner_fb fv_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4"

#FASTER
$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "fv_dynamics_mod.fv_dynamics (u v pt delp q w delz)/(u v pt delp q w delz)" fv_dynamics.F90 a2b_edge.F90 boundary.F90 dyn_core.F90 fv_arrays.F90 fv_arrays_nlm.F90 fv_cmp.F90 fv_control.F90 fv_fill.F90 fv_grid_utils.F90 fv_mapz.F90 fv_nesting.F90 fv_sg.F90 fv_tracer2d.F90 nh_core.F90 nh_utils.F90 sw_core.F90 tp_core.F90 MAPL/*.F90 tools/*.F90 fms/*.F90 -nocheckpoint "a2b_edge_mod.a2b_ord4_fb a2b_edge_mod.a2b_ord2 dyn_core_mod.dyn_core dyn_core_mod.pk3_halo dyn_core_mod.pln_halo dyn_core_mod.pe_halo dyn_core_mod.adv_pe dyn_core_mod.p_grad_c dyn_core_mod.nh_p_grad dyn_core_mod.split_p_grad dyn_core_mod.one_grad_p dyn_core_mod.grad1_p_update dyn_core_mod.mix_dp dyn_core_mod.geopk dyn_core_mod.del2_cubed dyn_core_mod.Rayleigh_fast fv_dynamics_mod.fv_dynamics fv_dynamics_mod.Rayleigh_Super fv_dynamics_mod.Rayleigh_Friction fv_dynamics_mod.compute_aam fv_dynamics_mod.geos_to_fv3 fv_dynamics_mod.fv3_to_geos fv_grid_utils_mod.cubed_to_latlon fv_grid_utils_mod.c2l_ord4 fv_grid_utils_mod.c2l_ord2 fv_mapz_mod.Lagrangian_to_Eulerian fv_mapz_mod.compute_total_energy fv_mapz_mod.pkez fv_mapz_mod.remap_z fv_mapz_mod.map_scalar_fb fv_mapz_mod.map1_ppm_fb fv_mapz_mod.mapn_tracer_fb fv_mapz_mod.map1_q2_fb fv_mapz_mod.remap_2d fv_mapz_mod.scalar_profile_fb fv_mapz_mod.cs_profile_fb fv_mapz_mod.cs_limiters fv_mapz_mod.ppm_profile fv_mapz_mod.ppm_limiters fv_mapz_mod.steepz fv_mapz_mod.rst_remap fv_mapz_mod.mappm fv_mapz_mod.moist_cv fv_mapz_mod.moist_cp fv_mapz_mod.map1_cubic fv_restart_mod.d2c_setup fv_tracer2d_mod.tracer_2d_1L fv_tracer2d_mod.tracer_2d fv_tracer2d_mod.tracer_2d_nested nh_core_mod.Riem_Solver3 nh_utils_mod.update_dz_c nh_utils_mod.update_dz_d nh_utils_mod.Riem_Solver_c nh_utils_mod.Riem_Solver3test nh_utils_mod.imp_diff_w nh_utils_mod.RIM_2D nh_utils_mod.SIM3_solver nh_utils_mod.SIM3p0_solver nh_utils_mod.SIM1_solver nh_utils_mod.SIM_solver nh_utils_mod.edge_scalar nh_utils_mod.edge_profile nh_utils_mod.nest_halo_nh sw_core_mod.c_sw sw_core_mod.d_sw  sw_core_mod.divergence_corner sw_core_mod.divergence_corner_nest sw_core_mod.d2a2c_vect sw_core_mod.fill3_4corners sw_core_mod.fill2_4corners sw_core_mod.fill_4corners sw_core_mod.xtp_u_fb sw_core_mod.ytp_v_fb sw_core_mod.compute_divergence_damping_fb sw_core_mod.smag_corner_fb tp_core_mod.mp_ghost_ew tp_core_mod.fv_tp_2d_fb tp_core_mod.copy_corners_fb tp_core_mod.xppm_fb tp_core_mod.yppm_fb tp_core_mod.deln_flux_fb a2b_edge_mod.extrap_corner_fb fv_grid_utils_mod.great_circle_dist sw_core_mod.edge_interpolate4"

cd adm/
rename mod_diff.f90 adm.F90 *mod_diff.f90

#Remove files we do not need
rm -f mpp* init_hydro_adm.F90 fv_arrays_adm.F90 fv_arrays_nlm_adm.F90 fms_io_adm.F90

#Rename routines that are duplicated for FWD/BWD versions
sed -i -- 's/_FB_FWD/_FWD/g' *0
sed -i -- 's/_FB_BWD/_BWD/g' *0
sed -i -- 's/_fb_fwd = /_fwd = /g' *0 #Get functions

#Delete lines with "!UNCOM" as these are TLM specific
sed -i '/\!UNCOM/d' *_adm.F90

#Get rid of push/pop of pointers

echo " "
echo "REPLACING CHECKPOINTING OF POINTERS WITH JUST RESETTING THE POINTER"
echo "-------------------------------------------------------------------"

rm -f pointers.txt

sed -i -- 's/USE ISO_C_BINDING/\!USE ISO_C_BINDING/g' *0
sed -i -- 's/USE ADMM_TAPENADE_INTERFACE/\!USE ADMM_TAPENADE_INTERFACE/g' *0
sed -i -- 's/TYPE(C_PTR) :: cptr/\!TYPE(C_PTR) :: cptr/g' *0
sed -i -- 's/INTEGER :: unknown_shape_in/\!INTEGER :: unknown_shape_in/g' *0
sed -i -- 's/CALL PUSHPOINTER/\!CALL PUSHPOINTER/g' *0
sed -i -- 's/CALL POPPOINTER/\!CALL POPPOINTER/g' *0

grep -i 'CALL C_F_POINTER' *0 | cut -d ',' -f 2- | cut -d ',' -f 1 | awk '\!seen[$0]++' > pointers.txt

foreach line ( "`cat pointers.txt`" )
 set pointto = `grep -i "=>" *0 | grep -i -m 1 "$line " | cut -d ' ' -f 7-`
 set replace = `echo "CALL C_F_POINTER(cptr, "$line", "`
 set with = `echo $line "=>" $pointto \!`
 echo "Replace $replace" 
 echo "with $with"
 echo " "
 sed -i -- "s/${replace}/${with}/g" *0
end
rm -f pointers.txt

echo " "


#Create list of subroutines that need to be copied later
grep -i "routine" *0 | grep -i "_adm(" > ReplaceThese.txt

#Make copies that wont use the checkpoint interface
mkdir origtap
cp *0 origtap/


#Switch to our checkpointing interface
#-------------------------------------

##Switch reals to custom Tapenade interface
sed -i -- 's/CALL PUSHREAL4ARRAY/CALL PUSHREALARRAY/g' *0
sed -i -- 's/CALL POPREAL4ARRAY/CALL POPREALARRAY/g' *0
sed -i -- 's/CALL PUSHREAL4(/CALL PUSHREALARRAY(/g' *0
sed -i -- 's/CALL POPREAL4(/CALL POPREALARRAY(/g' *0

sed -i -- 's/CALL PUSHREAL8ARRAY/CALL PUSHREALARRAY/g' *0
sed -i -- 's/CALL POPREAL8ARRAY/CALL POPREALARRAY/g' *0
sed -i -- 's/CALL PUSHREAL8(/CALL PUSHREALARRAY(/g' *0
sed -i -- 's/CALL POPREAL8(/CALL POPREALARRAY(/g' *0

#Integers
sed -i -- 's/CALL PUSHINTEGER4ARRAY/CALL PUSHINTEGER/g' *0
sed -i -- 's/CALL PUSHINTEGER4/CALL PUSHINTEGER/g' *0
sed -i -- 's/CALL POPINTEGER4ARRAY/CALL POPINTEGER/g' *0
sed -i -- 's/CALL POPINTEGER4/CALL POPINTEGER/g' *0

#Controls
sed -i -- 's/CALL PUSHCONTROL1B(/CALL PUSHCONTROL(1,/g' *0
sed -i -- 's/CALL PUSHCONTROL2B(/CALL PUSHCONTROL(2,/g' *0
sed -i -- 's/CALL PUSHCONTROL3B(/CALL PUSHCONTROL(3,/g' *0
sed -i -- 's/CALL PUSHCONTROL4B(/CALL PUSHCONTROL(4,/g' *0
sed -i -- 's/CALL PUSHCONTROL5B(/CALL PUSHCONTROL(5,/g' *0
sed -i -- 's/CALL POPCONTROL1B(/CALL POPCONTROL(1,/g' *0
sed -i -- 's/CALL POPCONTROL2B(/CALL POPCONTROL(2,/g' *0
sed -i -- 's/CALL POPCONTROL3B(/CALL POPCONTROL(3,/g' *0
sed -i -- 's/CALL POPCONTROL4B(/CALL POPCONTROL(4,/g' *0
sed -i -- 's/CALL POPCONTROL5B(/CALL POPCONTROL(5,/g' *0


#Switch to our checkpointing interface but _adm verion and only reals
#--------------------------------------------------------------------

cd origtap/

#Switch reals to custom Tapenade interface (integers and pop not important as always the same precision)
sed -i -- 's/CALL PUSHREAL4ARRAY/CALL PUSHREALARRAY_ADM/g' *0
sed -i -- 's/CALL POPREAL4ARRAY/CALL POPREALARRAY_ADM/g' *0
sed -i -- 's/CALL PUSHREAL4(/CALL PUSHREALARRAY_ADM(/g' *0
sed -i -- 's/CALL POPREAL4(/CALL POPREALARRAY_ADM(/g' *0

#Now bring back the _adm routines only, these must still use the Tapenade checkpointing
#--------------------------------------------------------------------------------------

cd ../

foreach line ( "`cat ReplaceThese.txt`" )

  set file_sen = `echo $line | cut -d ':' -f 1`
  set file_get = `echo origtap/$file_sen`

  set routine = `echo $line | cut -d ' ' -f 3 | cut -d '(' -f 1`
  set lsstr = `echo "SUBROUTINE" $routine"("`
  set lestr = `echo "END SUBROUTINE" $routine`

  set ls_sen = `grep -n "$lsstr" $file_sen | cut -d ':' -f 1`
  set le_sen = `grep -n "$lestr" $file_sen | cut -d ':' -f 1`
  set ls_get = `grep -n "$lsstr" $file_get | cut -d ':' -f 1`
  set le_get = `grep -n "$lestr" $file_get | cut -d ':' -f 1`

  set lsm1_sen = `expr $ls_sen - 1`
  set lep1_sen = `expr $le_sen + 1`

  awk "NR<=${lsm1_sen}" $file_sen > a.txt
  awk "NR>=${ls_get} && NR<=${le_get}" $file_get > b.txt
  awk "NR>=${lep1_sen}" $file_sen > c.txt

  rm -rf ${file_sen}
  cat a.txt b.txt c.txt > ${file_sen}

  rm -rf a.txt b.txt c.txt

end

rm -rf origtap
rm -f ReplaceThese.txt


#Rename _FB subroutines and rename calls
grep -i '_FB(' *0 | grep -i 'SUBROUTINE' > LinesToAdapt.txt
foreach line ( "`cat LinesToAdapt.txt`" )
 set file = `echo $line | cut -d ':' -f 1`
 set sub_fb = `echo $line | cut -d ' ' -f 3 | cut -d '(' -f 1`
 set sub_nl = `echo $sub_fb | rev | cut -c 4- | rev`
 set sub_fb_lc = `echo $sub_fb | tr 'A-Z' 'a-z'`
 set sub_nl_lc = `echo $sub_nl | tr 'A-Z' 'a-z'`
 set s1 = `echo SUBROUTINE ${sub_fb}\(`
 set s2 = `echo END SUBROUTINE ${sub_fb}`
 set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
 set le = `grep -in "$s2" $file | cut -d ':' -f 1`
 awk "NR<${ls}" $file > a.txt
 awk "NR>${le}" $file > b.txt
 cat a.txt b.txt > ${file}
 rm -f a.txt b.txt
 set r1 = `echo CALL ${sub_fb}\(`
 set r2 = `echo CALL ${sub_nl}\(`
 sed -i -- "s/${r1}/${r2}/g" *0
 sed -i -- "s/${sub_fb_lc}/${sub_nl_lc}/g" *0
end

rm -f LinesToAdapt.txt

#Rename _FB functions and rename calls
grep -i '_FB(' *0 | grep -i 'FUNCTION' > LinesToAdapt.txt
foreach line ( "`cat LinesToAdapt.txt`" )
 set file = `echo $line | cut -d ':' -f 1`
 set sub_fb = `echo $line | cut -d ' ' -f 4 | cut -d '(' -f 1`
 set sub_nl = `echo $sub_fb | rev | cut -c 4- | rev`
 set sub_fb_lc = `echo $sub_fb | tr 'A-Z' 'a-z'`
 set sub_nl_lc = `echo $sub_nl | tr 'A-Z' 'a-z'`
# exit()
 set s1 = `echo FUNCTION ${sub_fb}\(`
 set s2 = `echo END FUNCTION ${sub_fb}`
 set ls = `grep -in "$s1" $file | cut -d ':' -f 1`
 set le = `grep -in "$s2" $file | cut -d ':' -f 1`
# exit()
 awk "NR<${ls}" $file > a.txt
 awk "NR>${le}" $file > b.txt
 cat a.txt b.txt > ${file}
 rm -f a.txt b.txt
 set r1 = `echo = ${sub_fb}\(`
 set r2 = `echo = ${sub_nl}\(`
 sed -i -- "s/${r1}/${r2}/g" *0
 sed -i -- "s/${sub_fb_lc}/${sub_nl_lc}/g" *0
end

rm -f LinesToAdapt.txt



cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "Adjoint generation took $diff seconds"

#Done
