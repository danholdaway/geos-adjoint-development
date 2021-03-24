close all
clear
clc


cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_free_nocloud/
file = 'x0011dh_a.prog.eta.20130116_03z.nc4';

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_replay10_nocloud/
file = 'x0011dh_a.prog.eta.20130116_03z.nc4';

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

u_tlm = u_replay - u_free ; 
v_tlm = v_replay - v_free ; 
t_tlm = t_replay - t_free ; 
q_tlm = q_replay - q_free ; 
p_tlm = p_replay - p_free ; 
qi_tlm = qi_replay - qi_free ; 
ql_tlm = ql_replay - ql_free ; 
o3_tlm = o3_replay - o3_free ;

cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_replay10_nocloud/tlm_restart03/

file = 'fvpert.eta.nc4';
ncwrite(file,'u',u_tlm);
ncwrite(file,'v',v_tlm);
ncwrite(file,'tv',t_tlm);
ncwrite(file,'sphu',q_tlm);
ncwrite(file,'delp',p_tlm);
ncwrite(file,'qltot',ql_tlm);
ncwrite(file,'qitot',qi_tlm);
ncwrite(file,'ozone',o3_tlm);

% %Next
% 
% cd /discover/nobackup/drholdaw/tmp.22292/prog_replay10/
% file = 'x0011dh_a.prog.eta.20130116_01z.nc4';
% 
% u_replay = ncread(file,'u');
% v_replay = ncread(file,'v');
% t_replay = ncread(file,'tv');
% q_replay = ncread(file,'sphu');
% p_replay = ncread(file,'delp');
% qi_replay = ncread(file,'qitot');
% ql_replay = ncread(file,'qltot');
% o3_replay = ncread(file,'ozone');
% 
% u_tlm = u_replay - u_free ; 
% v_tlm = v_replay - v_free ; 
% t_tlm = t_replay - t_free ; 
% q_tlm = q_replay - q_free ; 
% p_tlm = p_replay - p_free ; 
% qi_tlm = qi_replay - qi_free ; 
% ql_tlm = ql_replay - ql_free ; 
% o3_tlm = o3_replay - o3_free ;
% 
% cd /discover/nobackup/drholdaw/tmp.22292/prog_replay10/tlm_restart
% 
% file = 'fvpert.eta.nc4';
% ncwrite(file,'u',u_tlm);
% ncwrite(file,'v',v_tlm);
% ncwrite(file,'tv',t_tlm);
% ncwrite(file,'sphu',q_tlm);
% ncwrite(file,'delp',p_tlm);
% ncwrite(file,'qltot',ql_tlm);
% ncwrite(file,'qitot',qi_tlm);
% ncwrite(file,'ozone',o3_tlm);
% 
% cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/BoundaryLayer/