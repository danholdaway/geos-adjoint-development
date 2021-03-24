close all
clear
clc

T = 255:0.001:285;

F2 = zeros(1,length(T));
DF2DT = zeros(1,length(T));
F2a = zeros(1,length(T));
F3 = zeros(1,length(T));

TE0=273.0;
TE1=263.0;
TE2=200.0;
RATE2 = 1.3;
RATE3 = 1.3;
JUMP1 = (RATE2-1.0) / ( ( TE0-TE1 )^0.333 ) ;

for i = 1:length(T);

    if ( T(i) >= TE0 )
        F2(i)    = 1.0;
        F2a(i)    = 1.0;
        F3(i)    = 1.0;
        DF2DT(i) = 0.0;
        DF2DTa(i) = 0.0;
    end

    if ( ( T(i) >= TE1 ) && ( T(i) < TE0 ) )
        F2(i)   = 1.0 + JUMP1 * (( TE0 - T(i) )^0.3333);
        F2a(i)   = 1.0 - (RATE3 - 1.0)/(TE0 - TE1) * (T(i) - TE0);
        F3(i)   = 1.0;
        DF2DT(i) =  -(0.3333*JUMP1)/(TE0-T(i))^0.6667;
        DF2DTa(i) = - (RATE3 - 1.0)/(TE0 - TE1);
    end

    if ( T(i) < TE1 ) 
        F2(i)   = RATE2 ;
        F2a(i)   = RATE2 ;
        F3(i)   = 1.0;
        DF2DT(i) = 0.0;
        DF2DTa(i) = 0.0;
    end


end


plot(T,F2)
hold on
plot(T,F2a,'r')

figure
plot(T,DF2DT)
hold on
plot(T,DF2DTa,'r')

% figure
% plot(T,F3)