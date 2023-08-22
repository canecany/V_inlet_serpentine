close all
clear all
clc
% Code from Chainan Sailabada and Juan Ordonez
% version 08/22/2023
% function Cold plate heat exchanger for the PEBBS + example

Tmax = 125 ;     % Design temperature(maximum) of the cold plate degC
Tin =30;       % Temperature inlet of the coolant degC
W =  0.5;       % m  PEBBS deimension width
L = 0.1953 ;       % m  PEBBS dimension length
qc = 1780 ;     % W (total Heat loss)
td = 0.25;       % pipe diameter inch 3/4
nt = 3 ;        % number of serpentine turn
t = 1.5;        % thickness of the cold plate in 

% heat function example
[dtfinal,vreq,massreq,nm,Tr,RSA,RE] = heatf(Tmax,Tin,W,L,t,qc,td,nt)
%massreq = mass req for 1 channel[kg/s]
%RSA = Thermal Resistance of the cold plate K/W

% pressure function example
if nm ~= 0
    Dt = td;
    [tot_mass_req,dP] = pressure_drop(vreq,L,td,nt,Tr,nm,Tin)
else
    dP = 1e6;
    print('Conditon is not satisfied')
end



