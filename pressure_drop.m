function [mass_inlet,dP_tot] = pressure_drop(Vin,L,Dt,nt,Tr,nm,Tin)

% input variable
% V = velocity inlet m/s getting from the function heat
% L = Length of the cold plate m
% W = width of the cold plate m
% Dt = Tube diameter m
% nt = number of turn
% Tr = turn radius m
% nm = number of manifold

% output variable
% dP = pressure drop of the system of one 1 pebbs 
Trm = Tr;
D = Dt*25.4e-3;

T0 = 273.15+Tin; % Temperature inlet
ele = 'Water';
P0 = 101350 ;   % pressure
[rho, mu, nu, cp_liquid, lambda,Pr,Cond_f] = fluid_properties_code(T0,P0,ele);


Da = 1/4*pi()*D^2;
V = Vin;
Re = rho*V*D/mu

% Pressure section for 1 serpentine 1 set
if Re < 2300
    f = 64/Re;
else
    f= 0.316/(Re^(1/4));
end
dP_major = f*((L*(nt+1))/D)*rho*(V^2)/2;
K_minor = 0.3*(nt)*1.25;
dP_minor = K_minor*rho*(V^2)/2;
dP_tot = (dP_major+dP_minor);
mass_inlet = rho*Da*V*nm;
end
