clear
close all
clc
mydir = pwd;

cd /home/drholdaw/LinearisedPhysics/Inputs/

pref

cd(mydir)

pp = p_ref;
pp_pert = pp + 10*rand(1,73);


minrhcrit = 0.93;
tempmaxrh = 0.94;
turnrhcrit = 750;
pi_0 = 4*atan(1);

a1 = minrhcrit + (tempmaxrh-minrhcrit)/19.*((atan((2.*(pp-turnrhcrit)/(1020.-turnrhcrit)-1.)*tan(20.*pi_0/21.-0.5*pi_0))+0.5*pi_0)*21./pi_0-1.);
a1_pert = minrhcrit + (tempmaxrh-minrhcrit)/19.*((atan((2.*(pp_pert-turnrhcrit)/(1020.-turnrhcrit)-1.)*tan(20.*pi_0/21.-0.5*pi_0))+0.5*pi_0)*21./pi_0-1.);


alpha = 1-a1;
alpha_pert = 1-a1_pert;

plot(alpha,pp)
set(gca,'YDir','reverse')

hold on
plot(alpha,pp_pert,'r')

figure
plot(alpha,1:73)
set(gca,'YDir','reverse')

QT = 1.5e-3;
qs = 1.55e-3;

DQ = 2*alpha(70)*qs; 

QMX = QT + 0.5*DQ;
QMN = QT - 0.5*DQ;

qs = 1.5e-3:1e-5:1.6e-3;

CF = (QT + (alpha(70) - 1).*qs)./(2*alpha(70)*qs);

% plot(qs,CF)