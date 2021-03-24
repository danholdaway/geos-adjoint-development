close all
clear
clc
mydir = pwd;

cd ../Inputs

pref

cd(mydir)

p = p_ref;
clear p_ref mydir


kappa = 0.7;

ph = 0.5*(p(1:end-1) + p(2:end));
ptothekappa1 = ph.^kappa;

delpk = p(2:end).^kappa - p(1:end-1).^kappa;
dellogp= log(p(2:end)) - log(p(1:end-1));
ptothekappa2 = delpk ./ (kappa .* dellogp);


plot(ptothekappa1,1:72)
hold on
plot(ptothekappa2,1:72,'r--')
set(gca,'YDir','reverse')
ylim([1 72])