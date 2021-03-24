close all
clear
clc

heating_vs_heatingpert_container_6hr

R = R1;
P = P1;

k = 10000;

J1_max = zeros(k,12);
J2_max = zeros(k,12);
J3_max = zeros(k,12);
J4_max = zeros(k,12);

J1_max(1:length(R1),1) = R1(:,3);
J2_max(1:length(R1),1) = R1(:,4);
J3_max(1:length(R1),1) = P1(:,3);
J4_max(1:length(R1),1) = P1(:,4);

J1_max(1:length(R2),2) = R2(:,3);
J2_max(1:length(R2),2) = R2(:,4);
J3_max(1:length(R2),2) = P2(:,3);
J4_max(1:length(R2),2) = P2(:,4);

J1_max(1:length(R3),3) = R3(:,3);
J2_max(1:length(R3),3) = R3(:,4);
J3_max(1:length(R3),3) = P3(:,3);
J4_max(1:length(R3),3) = P3(:,4);

J1_max(1:length(R4),4) = R4(:,3);
J2_max(1:length(R4),4) = R4(:,4);
J3_max(1:length(R4),4) = P4(:,3);
J4_max(1:length(R4),4) = P4(:,4);

J1_max(1:length(R5),5) = R5(:,3);
J2_max(1:length(R5),5) = R5(:,4);
J3_max(1:length(R5),5) = P5(:,3);
J4_max(1:length(R5),5) = P5(:,4);

J1_max(1:length(R6),6) = R6(:,3);
J2_max(1:length(R6),6) = R6(:,4);
J3_max(1:length(R6),6) = P6(:,3);
J4_max(1:length(R6),6) = P6(:,4);

J1_max(1:length(R7),7) = R7(:,3);
J2_max(1:length(R7),7) = R7(:,4);
J3_max(1:length(R7),7) = P7(:,3);
J4_max(1:length(R7),7) = P7(:,4);

J1_max(1:length(R8),8) = R8(:,3);
J2_max(1:length(R8),8) = R8(:,4);
J3_max(1:length(R8),8) = P8(:,3);
J4_max(1:length(R8),8) = P8(:,4);

J1_max(1:length(R9),9) = R9(:,3);
J2_max(1:length(R9),9) = R9(:,4);
J3_max(1:length(R9),9) = P9(:,3);
J4_max(1:length(R9),9) = P9(:,4);

J1_max(1:length(R10),10) = R10(:,3);
J2_max(1:length(R10),10) = R10(:,4);
J3_max(1:length(R10),10) = P10(:,3);
J4_max(1:length(R10),10) = P10(:,4);

J1_max(1:length(R11),11) = R11(:,3);
J2_max(1:length(R11),11) = R11(:,4);
J3_max(1:length(R11),11) = P11(:,3);
J4_max(1:length(R11),11) = P11(:,4);

J1_max(1:length(R12),12) = R12(:,3);
J2_max(1:length(R12),12) = R12(:,4);
J3_max(1:length(R12),12) = P12(:,3);
J4_max(1:length(R12),12) = P12(:,4);

count_J1 = zeros(1,12);
count_J2 = zeros(1,12);
count_J3 = zeros(1,12);
count_J4 = zeros(1,12);

for j = 1:12
    for i = 1:k

        if J4_max(i,j) > 10^(-3.5)

            count_J4(j) = count_J4(j) + 1;

        end

    end
end