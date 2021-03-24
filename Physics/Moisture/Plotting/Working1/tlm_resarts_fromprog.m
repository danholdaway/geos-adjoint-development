close all
clear
clc


cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free/
file = 'v000_C180.prog.eta.20140131_2200z.nc4';

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');
cfcn_free = ncread(file,'cfcn');

cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay/
file = 'v000_C180.prog.eta.20140131_2300z.nc4';

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');
cfcn_replay = ncread(file,'cfcn');

u_tlm = u_replay - u_free ; 
v_tlm = v_replay - v_free ; 
t_tlm = t_replay - t_free ; 
q_tlm = q_replay - q_free ; 
p_tlm = p_replay - p_free ; 
qi_tlm = qi_replay - qi_free ; 
ql_tlm = ql_replay - ql_free ; 
o3_tlm = o3_replay - o3_free ;
cfcn_tlm = cfcn_replay - cfcn_free ;

cd /discover/nobackup/drholdaw/btmp.25956/tlmrestarts_03z_nosinks/

file = 'fvpert.eta.nc4';
ncwrite(file,'u',u_tlm);
ncwrite(file,'v',v_tlm);
ncwrite(file,'tv',t_tlm);
ncwrite(file,'sphu',q_tlm);
ncwrite(file,'delp',p_tlm);
ncwrite(file,'qltot',ql_tlm);
ncwrite(file,'qitot',qi_tlm);
ncwrite(file,'ozone',o3_tlm);
ncwrite(file,'cfcn',cfcn_tlm);


file = 'fvpertX.eta.nc4';
ncwrite(file,'U',u_tlm);
ncwrite(file,'V',v_tlm);
ncwrite(file,'TV',t_tlm);
ncwrite(file,'QV',q_tlm);
ncwrite(file,'DP',p_tlm);
ncwrite(file,'QL',ql_tlm);
ncwrite(file,'QI',qi_tlm);
ncwrite(file,'O3',o3_tlm);
ncwrite(file,'CFCN',cfcn_tlm);

