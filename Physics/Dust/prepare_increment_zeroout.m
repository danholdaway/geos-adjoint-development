close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/wrk.dust/inc/

file_read = 'C360_513_snamma.inc.dust.eta.20060827_00z.nc4';

u = ncread(file_read,'u');

im = size(u,1);
jm = size(u,2);
lm = size(u,3);

u = zeros(im,jm,lm);
v = zeros(im,jm,lm);
t = zeros(im,jm,lm);
q = zeros(im,jm,lm);
p = zeros(im,jm,lm);
qi = zeros(im,jm,lm);
ql = zeros(im,jm,lm); 
o3 = zeros(im,jm,lm);
DU001 = zeros(im,jm,lm);
DU002 = zeros(im,jm,lm);
DU003 = zeros(im,jm,lm);
DU004 = zeros(im,jm,lm);
DU005 = zeros(im,jm,lm);

file_write = 'C360_513_snamma.inc.dust.eta.20060827_00z.nc4';

ncwrite(file_write,'u',u);
ncwrite(file_write,'v',v);
ncwrite(file_write,'tv',t);
ncwrite(file_write,'sphu',q);
ncwrite(file_write,'delp',p);
ncwrite(file_write,'qltot',ql);
ncwrite(file_write,'qitot',qi);
ncwrite(file_write,'ozone',o3);
ncwrite(file_write,'DU001',DU001);
ncwrite(file_write,'DU002',DU002);
ncwrite(file_write,'DU003',DU003);
ncwrite(file_write,'DU004',DU004);
ncwrite(file_write,'DU005',DU005);
