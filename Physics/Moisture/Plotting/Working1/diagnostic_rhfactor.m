close all
clear
clc

LM = 72

L = 1:1:LM;
RhFtop = 1.10

RhF = (-(RhFtop-1)/(LM-1))*L + 1 + ((RhFtop-1)/(LM-1))*LM;

plot(RhF,L)

set(gca,'YDir','reverse')
ylim([1 72])

a = 4263.521+3043.788+2369.931+2113.617+1981.326+1918.161+1896.582+1895.765+1901.893+1901.030+1892.357+1890.435+1915.948+1949.788+1968.891+1935.187+1847.901+1759.155+1686.419+1622.992+1560.990+1491.463+1433.044+1382.247+1336.741+1298.246+1256.677+1225.371+1196.935+1168.765+1141.295+1115.706+1092.122+1074.046+1051.416+1029.569+1008.620+999.7606+1008.404+1016.308+1024.119+1015.889+1056.928+1090.771+1101.740+755.8752+699.0190+648.5204+605.5278+568.0486+536.5005+508.9009+482.1696+457.7308+292.8661+284.0783+275.6737+267.5493+259.3953+202.1328+148.3225+145.2696+143.1035+141.1251+139.0071+137.4046+134.8894+132.9001+131.3135+129.8465+128.4180+126.1835

dzet = [4263.521 3043.788 2369.931 2113.617 1981.326 1918.161 1896.582 1895.765 1901.893 1901.030 1892.357 1890.435 1915.948 1949.788 1968.891 1935.187 1847.901 1759.155 1686.419 1622.992 1560.990 1491.463 1433.044 1382.247 1336.741 1298.246 1256.677 1225.371 1196.935 1168.765 1141.295 1115.706 1092.122 1074.046 1051.416 1029.569 1008.620 999.7606 1008.404 1016.308 1024.119 1015.889 1056.928 1090.771 1101.740 755.8752 699.0190 648.5204 605.5278 568.0486 536.5005 508.9009 482.1696 457.7308 292.8661 284.0783 275.6737 267.5493 259.3953 202.1328 148.3225 145.2696 143.1035 141.1251 139.0071 137.4046 134.8894 132.9001 131.3135 129.8465 128.4180 126.1835];

z_ref = zeros(1,LM+1);

z_ref(LM+1) = 0;
for i = LM:-1:1;
    z_ref(i) = z_ref(i+1) + dzet(i);
end
zrefh = 0.5*(z_ref(1:end-1)+z_ref(2:end));

cd ../../Inputs/
pref
prefh = 0.5*(p_ref(1:end-1)+p_ref(2:end));
cd ../DAS_plotting/Cloud' Fraction'/


figure
plot(dzet,L);
set(gca,'YDir','reverse')
ylim([1 72])


figure
plot(dzet,prefh,'-x');
set(gca,'YDir','reverse')

figure
plot(dzet,zrefh,'-x');
% set(gca,'YDir','reverse')