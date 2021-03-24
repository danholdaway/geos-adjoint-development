close all
clear
clc

cd ../Inputs/

pref

cd ../DAS_plotting/

pe = p_ref/0.01;

kappa = 5/7;

p0 = 1.0;
% p0 = 100000.0;

ph = 0.5*(pe(1:end-1) + pe(2:end));
con1 = (p0./ph).^kappa;



pk = pe.^kappa;

pke_local = log(pe);

bx_local = pke_local(2:end)-pke_local(1:end-1);

pke_local = pe.^kappa;

cx_local = pke_local(2:end)-pke_local(1:end-1);

bx_local = 1./(kappa*bx_local);

pkz = bx_local.*cx_local;


