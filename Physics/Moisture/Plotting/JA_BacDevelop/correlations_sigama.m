close all
clear
clc
mydir = pwd;


cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/Correlations/

cfnd_vs_sigmaqt1
qcnd_vs_sigmaqt1

cd(mydir)

figure

subplot(2,1,1)
scatter(cfnd,sigma1_c,'xk')

subplot(2,1,2)
scatter(qcnd,sigma1_q,'xk')