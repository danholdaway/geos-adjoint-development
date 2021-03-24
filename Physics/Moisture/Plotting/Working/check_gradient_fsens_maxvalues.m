close all
clear
clc

%ENTER DAY AND STARTING HOUR
d2 = '20130102';
H = '21';

daten = datenum(d2, 'yyyymmdd');
daten = daten-1;

u_max = zeros(1,31);
v_max = zeros(1,31);
t_max = zeros(1,31);
q_max = zeros(1,31);
p_max = zeros(1,31);
qi_max = zeros(1,31);
ql_max = zeros(1,31);
o3_max = zeros(1,31);

jmin = 1;
jmax = 31;

% cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/filt5_rmd5_adjd
cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/SBAB/FirstTest/

for j = jmin:jmax
   
    daten = daten + 1;

    d1 = datestr(daten-1, 'yyyymmdd');
    d2 = datestr(daten  , 'yyyymmdd');
    d3 = datestr(daten+1, 'yyyymmdd');

    disp([num2str(j) ' ' d2])
            
    file = ['x0011dh_a.fsens_tce.eta.',d1,'_',H,'z+',d3,'_00z-',d2,'_00z.nc4'];

    lon = ncread(file,'lon');
    lat = ncread(file,'lat');
    u = ncread(file,'u');
    v = ncread(file,'v');
    t = ncread(file,'tv');
    q = ncread(file,'sphu');
    p = ncread(file,'delp');
    qi = ncread(file,'qitot');
    ql = ncread(file,'qltot');
    o3 = ncread(file,'ozone');

    u_max(j) = max(abs(u(:)));
    v_max(j) = max(abs(v(:)));
    t_max(j) = max(abs(t(:)));
    q_max(j) = max(abs(q(:)));
    p_max(j) = max(abs(p(:)));
    qi_max(j) = max(abs(qi(:)));
    ql_max(j) = max(abs(ql(:)));
    o3_max(j) = max(abs(o3(:)));

end

disp('u')
disp([(jmin:jmax)' u_max(jmin:jmax)'])
disp('v')
disp([(jmin:jmax)' v_max(jmin:jmax)'])
disp('t')
disp([(jmin:jmax)' t_max(jmin:jmax)'])
disp('q')
disp([(jmin:jmax)' q_max(jmin:jmax)'])
disp('p')
disp([(jmin:jmax)' p_max(jmin:jmax)'])
disp('qi')
disp([(jmin:jmax)' qi_max(jmin:jmax)'])
disp('ql')
disp([(jmin:jmax)' ql_max(jmin:jmax)'])
disp('o3')
disp([(jmin:jmax)' o3_max(jmin:jmax)'])



cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/
