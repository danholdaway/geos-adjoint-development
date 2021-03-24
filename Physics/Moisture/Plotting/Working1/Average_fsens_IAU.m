close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/tmp.20094/sens.20130104.000000

im = 576;
jm = 361;
lm = 72;

u = zeros(im,jm,lm);
v = zeros(im,jm,lm);
t = zeros(im,jm,lm);
q = zeros(im,jm,lm);
p = zeros(im,jm,lm);

for k = 1001:1024
    
    
    file = ['fsens.iau.',num2str(k),'.nc4'];

    disp(file)
    
    u = u + ncread(file,'u');
    v = v + ncread(file,'v');
    t = t + ncread(file,'tv');
    q = q + ncread(file,'sphu');
    p = p + ncread(file,'delp');
    
end
    
u = u / 24;
v = v / 24;
t = t / 24;
q = q / 24;
p = p / 24;

file_write = 'fsens.iau.average.nc4';

ncwrite(file_write,'u',u);
ncwrite(file_write,'v',v);
ncwrite(file_write,'t',t);
ncwrite(file_write,'q',q);
ncwrite(file_write,'p',p);


cd(mydir)