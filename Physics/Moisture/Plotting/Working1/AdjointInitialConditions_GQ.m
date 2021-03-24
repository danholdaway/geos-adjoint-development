close all
clear
clc

% Code to read the Jgradf.eta.nc4 from two directories and average them.
% This provides the initital conditions for the adjoint when running 
% with Gaussian Quadrature.

% Note that Jgradf.eta.nc4 needs to be copied to JgradfX.eta.nc4,
% variables need to be renamed.

% Copy Jgradf.eta.nc4 to JgradfGQ.eta.nc4 in the first directory before
% running

% Names of the files containing adjoint restarts
file_read  = 'Jgradf.eta.nc4';
file_write = 'JgradfGQ.eta.nc4';

% Move to directory where FIRST forecast creates adjoint initial conditions
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/

u_a     = ncread(file_read,'u');
v_a     = ncread(file_read,'v');
tv_a    = ncread(file_read,'tv');
sphu_a  = ncread(file_read,'sphu');
delp_a  = ncread(file_read,'delp');
qitot_a = ncread(file_read,'qitot');
qltot_a = ncread(file_read,'qltot');
ozone_a = ncread(file_read,'ozone');

% Move to directory where SECOND forecast creates adjoint initial conditions
cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000/

u_b     = ncread(file_read,'u');
v_b     = ncread(file_read,'v');
tv_b    = ncread(file_read,'tv');
sphu_b  = ncread(file_read,'sphu');
delp_b  = ncread(file_read,'delp');
qitot_b = ncread(file_read,'qitot');
qltot_b = ncread(file_read,'qltot');
ozone_b = ncread(file_read,'ozone');

% Average the fields
u_c     = 0.5*(u_a     + u_b    );
v_c     = 0.5*(v_a     + v_b    );
tv_c    = 0.5*(tv_a    + tv_b   );
sphu_c  = 0.5*(sphu_a  + sphu_b );
delp_c  = 0.5*(delp_a  + delp_b );
qitot_c = 0.5*(qitot_a + qitot_b);
qltot_c = 0.5*(qltot_a + qltot_b);
ozone_c = 0.5*(ozone_a + ozone_b);

% Move back to FIRST forecast location
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/

% Write averaged fields.
ncwrite(file_write,'u',u_c);
ncwrite(file_write,'v',v_c);
ncwrite(file_write,'tv',tv_c);
ncwrite(file_write,'sphu',sphu_c);
ncwrite(file_write,'delp',delp_c);
ncwrite(file_write,'qitot',qitot_c);
ncwrite(file_write,'qltot',qltot_c);
ncwrite(file_write,'ozone',ozone_c);



