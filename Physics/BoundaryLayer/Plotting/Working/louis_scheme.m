close all
clc
clear

%Choose model level to plot.
plot_level = 50;

%Choose whether to match the scale.
match_scale = 1;

cd /discover/nobackup/drholdaw/tmp.22292/BLexp/control/
file = 'x0011dh_a.prog.eta.20130116_06z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');

u = ncread(file,'u');
v = ncread(file,'v');
th = ncread(file,'tv');
q = ncread(file,'sphu');
dp = ncread(file,'delp');
qi = ncread(file,'qitot');
ql = ncread(file,'qltot');
o3 = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/BoundaryLayer/

[im,jm,lm] = size(u);

P = zeros(im,jm,lm+1);

P(:,:,1) = 1;

for l = 1:lm
    P(:,:,l+1) = dp(:,:,l) + P(:,:,l);
end
% P = P*0.01;
Pf = 0.5*(P(:,:,1:end-1) + P(:,:,2:end));

CONS_P00 = 100000.0;
CONS_RUNIV = 8314.3;
CONS_KAPPA = 2.0/7.0;
CONS_AIRMW = 28.97;
CONS_H2OMW = 18.01;
CONS_GRAV = 9.80;
CONS_KARMAN = 0.40;
CONS_TICE   = 273.16;
CONS_ALHL = 2.4665E6;
CONS_ALHF = 3.3370E5;
CONS_ALHS = CONS_ALHL + CONS_ALHF;
CONS_RGAS = CONS_RUNIV/CONS_AIRMW;
CONS_CP = CONS_RGAS/CONS_KAPPA;
CONS_VIREPS = CONS_AIRMW/CONS_H2OMW-1.0;

%Compute the edge heights using Arakawa-Suarez hydrostatic equation
Z = zeros(im,jm,lm+1);

Z(:,:,lm+1) = 0.0;

PKE_atLp1 = zeros(im,jm);
PKE_atL = zeros(im,jm);

for L = lm:-1:1
   PKE_atLp1  = (P(:,:,L+1)/CONS_P00).^CONS_KAPPA;
   PKE_atL    = (P(:,:,L  )/CONS_P00).^CONS_KAPPA;

   Z(:,:,L) = Z(:,:,L+1) + (CONS_CP/CONS_GRAV).*th(:,:,L).*(PKE_atLp1-PKE_atL);
end

Zf = 0.5*(Z(:,:,1:end-1) + Z(:,:,2:end));

%Exner pressures
pi  = (P /CONS_P00).^(CONS_RGAS/CONS_CP);
pif = (Pf/CONS_P00).^(CONS_RGAS/CONS_CP);

t = pif.*th;
tv  = t .* ( 1.0 + CONS_VIREPS * q - ql - qi );

%Virtual potential temperature with 1-2-1 smooth of bottom 5 levels
thv = tv.*(th./t);
thv(:,:,lm  ) = thv(:,:,lm-1)*0.25 + thv(:,:,lm  )*0.75;
thv(:,:,lm-1) = thv(:,:,lm-2)*0.25 + thv(:,:,lm-1)*0.50 + thv(:,:,lm  )*0.25 ;
thv(:,:,lm-2) = thv(:,:,lm-3)*0.25 + thv(:,:,lm-2)*0.50 + thv(:,:,lm-1)*0.25 ;
thv(:,:,lm-3) = thv(:,:,lm-4)*0.25 + thv(:,:,lm-3)*0.50 + thv(:,:,lm-2)*0.25 ;
thv(:,:,lm-4) = thv(:,:,lm-5)*0.25 + thv(:,:,lm-4)*0.50 + thv(:,:,lm-3)*0.25 ;
thv(:,:,lm-5) = thv(:,:,lm-6)*0.25 + thv(:,:,lm-5)*0.50 + thv(:,:,lm-4)*0.25 ;




