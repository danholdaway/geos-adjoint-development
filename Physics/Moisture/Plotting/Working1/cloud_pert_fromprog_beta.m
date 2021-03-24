close all
clear
clc


cd /discover/nobackup/drholdaw/btmp.25956/hord_tr/hord_tr12/prog_free/
file = 'v000_C180.prog.eta.20140201_0000z.nc4';

qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

qi_tlm = 1.1*qi_free - qi_free;
ql_tlm = 1.1*ql_free - ql_free;

im = size(qi_free,1);
jm = size(qi_free,2);
lm = size(qi_free,3);

u_tlm = zeros(im,jm,lm);
v_tlm = zeros(im,jm,lm);
t_tlm = zeros(im,jm,lm);
q_tlm = zeros(im,jm,lm);
p_tlm = zeros(im,jm,lm);
o3_tlm = zeros(im,jm,lm);


cd /discover/nobackup/drholdaw/btmp.25956/hord_tr/hord_tr12/tlmrestarts_00z_cloudonly_b10/
file = 'fvpert.eta.nc4';
ncwrite(file,'u',u_tlm);
ncwrite(file,'v',v_tlm);
ncwrite(file,'tv',t_tlm);
ncwrite(file,'sphu',q_tlm);
ncwrite(file,'delp',p_tlm);
ncwrite(file,'qltot',ql_tlm);
ncwrite(file,'qitot',qi_tlm);
ncwrite(file,'ozone',o3_tlm);



