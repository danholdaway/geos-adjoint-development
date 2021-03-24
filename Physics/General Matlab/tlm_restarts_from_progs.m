close all
clear
clc


cd /discover/nobackup/drholdaw/wrk.sandy/prog/prog_free/
file = 'sens_sandy.prog.eta.20121029_03z.nc4';

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');
% DU001_free = ncread(file,'DU001');
% DU002_free = ncread(file,'DU002');
% DU003_free = ncread(file,'DU003');
% DU004_free = ncread(file,'DU004');
% DU005_free = ncread(file,'DU005');

cd /discover/nobackup/drholdaw/wrk.sandy/prog/prog_replay/
file = 'sens_sandy.prog.eta.20121029_03z.nc4';

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');
% DU001_replay = ncread(file,'DU001');
% DU002_replay = ncread(file,'DU002');
% DU003_replay = ncread(file,'DU003');
% DU004_replay = ncread(file,'DU004');
% DU005_replay = ncread(file,'DU005');

u_tlm = u_replay - u_free ; 
v_tlm = v_replay - v_free ; 
t_tlm = t_replay - t_free ; 
q_tlm = q_replay - q_free ; 
p_tlm = p_replay - p_free ; 
qi_tlm = qi_replay - qi_free ; 
ql_tlm = ql_replay - ql_free ; 
o3_tlm = o3_replay - o3_free ;
% DU001_tlm = DU001_replay - DU001_free ;
% DU002_tlm = DU002_replay - DU002_free ;
% DU003_tlm = DU003_replay - DU003_free ;
% DU004_tlm = DU004_replay - DU004_free ;
% DU005_tlm = DU005_replay - DU005_free ;

cd /discover/nobackup/drholdaw/wrk.sandy/tlm_restarts_20121029_030000/

file = 'fvpert.eta.nc4';
ncwrite(file,'u',u_tlm);
ncwrite(file,'v',v_tlm);
ncwrite(file,'tv',t_tlm);
ncwrite(file,'sphu',q_tlm);
ncwrite(file,'delp',p_tlm);
ncwrite(file,'qltot',ql_tlm);
ncwrite(file,'qitot',qi_tlm);
ncwrite(file,'ozone',o3_tlm);
% ncwrite(file,'DU001',DU001_tlm);
% ncwrite(file,'DU002',DU002_tlm);
% ncwrite(file,'DU003',DU003_tlm);
% ncwrite(file,'DU004',DU004_tlm);
% ncwrite(file,'DU005',DU005_tlm);

file = 'fvpertX.eta.nc4';
ncwrite(file,'U',u_tlm);
ncwrite(file,'V',v_tlm);
ncwrite(file,'TV',t_tlm);
ncwrite(file,'QV',q_tlm);
ncwrite(file,'DP',p_tlm);
ncwrite(file,'QL',ql_tlm);
ncwrite(file,'QI',qi_tlm);
ncwrite(file,'O3',o3_tlm);
% ncwrite(file,'DU001',DU001_tlm);
% ncwrite(file,'DU002',DU002_tlm);
% ncwrite(file,'DU003',DU003_tlm);
% ncwrite(file,'DU004',DU004_tlm);
% ncwrite(file,'DU005',DU005_tlm);
