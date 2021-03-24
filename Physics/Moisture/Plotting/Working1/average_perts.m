close all
clear
clc


cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = 'v000_C180.fvpert.eta.20140202_0300z_ST03z_P22000_S1.nc4';

u_free = ncread(file,'U');
v_free = ncread(file,'V');
t_free = ncread(file,'TV');
q_free = ncread(file,'QV');
p_free = ncread(file,'DP');
qi_free = ncread(file,'QI');
ql_free = ncread(file,'QL');
o3_free = ncread(file,'O3');

cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = 'v000_C180.fvpert.eta.20140202_0300z_ST03z_P22000_S1_TRAJ2.nc4';

u_replay = ncread(file,'U');
v_replay = ncread(file,'V');
t_replay = ncread(file,'TV');
q_replay = ncread(file,'QV');
p_replay = ncread(file,'DP');
qi_replay = ncread(file,'QI');
ql_replay = ncread(file,'QL');
o3_replay = ncread(file,'O3');

cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = 'v000_C180.fvpert.eta.20140202_0300z_ST03z_P22000_S1_MEANTRAJ.nc4';

u_tlm = 0.5*(u_replay + u_free) ; 
v_tlm = 0.5*(v_replay + v_free) ; 
t_tlm = 0.5*(t_replay + t_free) ; 
q_tlm = 0.5*(q_replay + q_free) ; 
p_tlm = 0.5*(p_replay + p_free) ; 
qi_tlm = 0.5*(qi_replay + qi_free) ; 
ql_tlm = 0.5*(ql_replay + ql_free) ; 
o3_tlm = 0.5*(o3_replay + o3_free) ;

ncwrite(file,'U',u_tlm);
ncwrite(file,'V',v_tlm);
ncwrite(file,'TV',t_tlm);
ncwrite(file,'QV',q_tlm);
ncwrite(file,'DP',p_tlm);
ncwrite(file,'QL',ql_tlm);
ncwrite(file,'QI',qi_tlm);
ncwrite(file,'O3',o3_tlm);


