close all
clear all
clc
% Code from Chainan Sailabada and Juan Ordonez
% version 06/08/2023
% function Cold plate heat exchanger for the PEBBS + example
n = 100;
for i =1:n
Tmax = 125 ;     % Design temperature(maximum) of the cold plate degC
Tin =linspace(15,90,n);       % Temperature inlet of the coolant degC
W =  0.5;       % m  PEBBS deimension width
L = 0.1953 ;       % m  PEBBS dimension length
qc = 1780 ;     % W (total Heat loss)
td = 0.25;       % pipe diameter inch 3/4
nt = 3 ;        % number of serpentine turn
t = 1.5;        % thickness of the cold plate in 

% heat function example
[dtfinal(i),vreq(i),massreq,nm,Tr,RSA,RE] = heatf(Tmax,Tin(i),W,L,t,qc,td,nt)
%massreq = mass req for 1 channel[kg/s]
%RSA = Thermal Resistance of the cold plate K/W

% pressure function example
if nm ~= 0
    Dt = td;
    [tot_mass_req,dP(i)] = pressure_drop(vreq(i),L,td,nt,Tr,nm,Tin(i))
else
    dP = 1e6
    print('Conditon is not satisfied')
end
end
% figure(1)
% subplot(2,2,1)
% plot(Tin,dP)
% xlabel('Tin C')
% ylabel('Pressure drop Pa')
% title('Temperature inlet and Pressure drop ')
% 
% subplot(2,2,2)
% plot(Tin,dtfinal)
% xlabel('Tin C')
% ylabel('Cold plate surface temperature C')
% title('Temperature inlet and surface temperature')
% 
% subplot(2,2,3)
% plot(Tin,vreq)
% xlabel('Tin C')
% ylabel('Velocity require m/s')
% title('Temperature inlet and Velocity require')











