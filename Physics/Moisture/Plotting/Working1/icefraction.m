close all
clear
clc

TEMP = 220:0.1:280;
ICEFRCT = zeros(1,length(TEMP));
ICEFRCT1 = zeros(1,length(TEMP));


T_ICE_ALL = 273.16 - 40;
T_ICE_MAX = 273.16;
ICEFRPWR = round( 4.0 + .001 );

for i = 1:length(TEMP);

    if ( TEMP(i) <= T_ICE_ALL )
        ICEFRCT(i) = 1.0;
    elseif ( (TEMP(i) > T_ICE_ALL) && (TEMP(i) <= T_ICE_MAX) )
        ICEFRCT(i) = 1.00 -  ( TEMP(i) - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) ;
    end
  
    if ( TEMP(i) <= T_ICE_ALL )
        ICEFRCT1(i) = 1.0;
    elseif ( (TEMP(i) > T_ICE_ALL) && (TEMP(i) <= T_ICE_MAX) )
        ICEFRCT1(i) =  1.00 -  ( TEMP(i) - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) ;
    end

    ICEFRCT(i) = min(ICEFRCT(i),1.0);
    ICEFRCT(i) = max(ICEFRCT(i),0.0);
        
end

ICEFRCT = ICEFRCT.^ICEFRPWR;
ICEFRCT1 = ICEFRCT1.^2;

plot(TEMP,ICEFRCT)
hold on
plot(TEMP,ICEFRCT1,'r')