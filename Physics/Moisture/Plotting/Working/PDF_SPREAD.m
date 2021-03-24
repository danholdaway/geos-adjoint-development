close all
clear
clc

pp = 1:0.1:1000;

a1 = zeros(1,length(pp));
alpha = zeros(1,length(pp));

maxrhcrit = 1.0;
minrhcrit = 0.93; %Resolution Dependent
maxrhcritland = minrhcrit + 0.01;
turnrhcrit = 750.0;
pi_0 = 4.*atan(1.);

for j = 1:2

    if j == 1
        frland = 1.0;
    else
        frland = 0.01;
    end

    tempmaxrh = maxrhcrit;

    if (frland > 0.05)
        tempmaxrh = maxrhcritland;
    end

    for i = 1:length(pp)

        a1(j,i) = 1.0;

        if (pp(i) <= turnrhcrit)

        a1(j,i) = minrhcrit;

        else

        a1(j,i) = minrhcrit + (tempmaxrh-minrhcrit)/(19.)*((atan( (2.*(pp(i)- turnrhcrit)/(1020.-turnrhcrit)-1.)*tan(20.*pi_0/21.-0.5*pi_0) ) + 0.5*pi_0) * 21./pi_0 - 1.);

        end

        a1(j,i) = min(a1(j,i),1.);

        alpha(j,i) = 1. - a1(j,i);

        alpha(j,i) = min( alpha(j,i) , 0.25 ); 

    end

end
    
figure
plot(a1(1,:),pp)
hold on
plot(a1(2,:),pp,'r')

figure
plot(alpha(1,:),pp)
hold on
plot(alpha(2,:),pp,'r')