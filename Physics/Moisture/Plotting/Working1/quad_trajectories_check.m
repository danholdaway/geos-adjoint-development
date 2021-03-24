close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/traj_halfdeg/
file = 'x0011dh_a.traj.lcv.20130116_0000z.nc4';

U_a = ncread(file,'U');
V_a = ncread(file,'V');
PT_a = ncread(file,'PT');
QV_a = ncread(file,'QV');
DP_a = ncread(file,'DP');
QI_a = ncread(file,'QI');
QL_a = ncread(file,'QL');
CFLS_a = ncread(file,'CFLS');
CFCN_a = ncread(file,'CFCN');
QILS_a = ncread(file,'QILS');
QICN_a = ncread(file,'QICN');
QLLS_a = ncread(file,'QLLS');
QLCN_a = ncread(file,'QLCN');
O3_a = ncread(file,'O3');
KCBL_a = ncread(file,'KCBL');
TS_a = ncread(file,'TS');

cd /discover/nobackup/drholdaw/tmp.22292/traj_quad/repl/

U_b = ncread(file,'U');
V_b = ncread(file,'V');
PT_b = ncread(file,'PT');
QV_b = ncread(file,'QV');
DP_b = ncread(file,'DP');
QI_b = ncread(file,'QI');
QL_b = ncread(file,'QL');
CFLS_b = ncread(file,'CFLS');
CFCN_b = ncread(file,'CFCN');
QILS_b = ncread(file,'QILS');
QICN_b = ncread(file,'QICN');
QLLS_b = ncread(file,'QLLS');
QLCN_b = ncread(file,'QLCN');
O3_b = ncread(file,'O3');
KCBL_b = ncread(file,'KCBL');
TS_b = ncread(file,'TS');

cd /discover/nobackup/drholdaw/tmp.22292/traj_quad/
U_c = ncread(file,'U');
V_c = ncread(file,'V');
PT_c = ncread(file,'PT');
QV_c = ncread(file,'QV');
DP_c = ncread(file,'DP');
QI_c = ncread(file,'QI');
QL_c = ncread(file,'QL');
CFLS_c = ncread(file,'CFLS');
CFCN_c = ncread(file,'CFCN');
QILS_c = ncread(file,'QILS');
QICN_c = ncread(file,'QICN');
QLLS_c = ncread(file,'QLLS');
QLCN_c = ncread(file,'QLCN');
O3_c = ncread(file,'O3');
KCBL_c = ncread(file,'KCBL');
TS_c = ncread(file,'TS');





