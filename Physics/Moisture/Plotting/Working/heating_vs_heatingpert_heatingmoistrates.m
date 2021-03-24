close all
clc
clear


cd /discover/nobackup/drholdaw/tmp.9748/sens.20130112.000000/

H_M_J


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


ivec = ivec - 1;

t = size(A,1);



for i = 1:t
    
   
    H1 = A(i,1:ivec(i),1);
    M1 = A(i,1:ivec(i),2);
    
    
    [H,ind] = sort(H1,'ascend');
    MbyH = M1(ind);
    [M,ind] = sort(M1,'ascend');
    HbyM = H1(ind);
    
    
    H_95(i) = H(round(95*length(H)/100));
    H_96(i) = H(round(96*length(H)/100));
    H_97(i) = H(round(97*length(H)/100));
    H_98(i) = H(round(98*length(H)/100));
    
    M_95(i) = M(round(95*length(M)/100));
    M_96(i) = M(round(96*length(M)/100));
    M_97(i) = M(round(97*length(M)/100));
    M_98(i) = M(round(98*length(M)/100));
    
end
   

semilogy(H)
hold on
semilogy(MbyH,'r')
plot(1:length(H),H_95(i)*ones(1,length(H)),'r')
plot(1:length(H),H_96(i)*ones(1,length(H)),'g')
plot(1:length(H),H_97(i)*ones(1,length(H)),'b')
plot(1:length(H),H_98(i)*ones(1,length(H)),'m')
xlim([1 length(H)])

figure

semilogy(M)
hold on
semilogy(HbyM,'r')
plot(1:length(M),M_95(i)*ones(1,length(M)),'r')
plot(1:length(M),M_96(i)*ones(1,length(M)),'g')
plot(1:length(M),M_97(i)*ones(1,length(M)),'b')
plot(1:length(M),M_98(i)*ones(1,length(M)),'m')
xlim([1 length(M)])


