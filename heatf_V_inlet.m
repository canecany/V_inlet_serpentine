function [dtfinal,massreq,nm,Tr,RSA,Re] = heatf_V_inlet(vin,Tin,W,L,tc,qc,td,nt)
%input variable
% vin = velocity inlet m/s
% Tin = temperature inlet degC
% W = width of PEEBS m
% L = Length of pebbs m
% tc = thickness oc cool plate %inch
% qc = heat form Pebbs W
% td = tube diameter %inch
% nt = number of turn inthe cold plate

%output variable
% dt final = temperature average of plate
% massreq = total mass flowreate require for the system
% nm = number of manifold system
% Tr = turn radius use for the pressure drop alculation
% RSA = thermal resistance of cold plate

T0 = 273+Tin;           % Temperature inlet
P0 = 101325;            % Pa
ele ="Water";           % define fluid element(such as "Water", "Helium","Nitrogen", any)
[rho, mu, nu, cp_liquid, lambda,Pr,Cond_f] = fluid_properties(T0,P0,ele);

% Geometry define
inM = 0.0254; % inch to meter unit converter
t = tc*inM;
k_s = 15 ;% thermal conductivity property of aluminum
%tube size 

Dm = td*inM;    % diameter m
wc = t*1.5;
Tr = wc/2;
%width chanael
nc = floor(W/wc);
nm = floor(nc/nt); %number of manifold

%Heat transfer area 
At   = pi()*Dm*L*nm*nc;        % area of heat transfer in one section
Acss = pi()*(Dm^2)/4;          % area coross section of tube
q_s = qc/nm ;                   % heat in one maniflow set
HTA = pi()*Dm*L*(nt+1);            % heat transfer area
tsub = 1/2*(1/2*(t-Dm)+1/4*Dm);
Low_Re = 2500;              %Heat exchanger need to be turbulentOke

Re = rho*vin*Dm/mu;                                             %correct
h = heat_transfer_coefficent(Pr,Re,Cond_f,Dm);                     %correct
massf = rho*Acss*vin;

Rconv = 1/(h*HTA);
Rthppbulk = (L-Tr*(nt/2))*(nt+1)*(2*pi*Dm)/(massf*cp_liquid);    
Rcond1 = ((pi*Dm)/(pi*k_s)) * log (1/(sin(0.25*pi*Dm/(pi*Dm))));
S = 2*pi*(L-Tr*(nt/2))*(nt+1)/(log(1.08*t/Dm));
Rcond2 = 1/(S*k_s);
Rthppcnd = tsub/k_s;
%Rtt = Rconv+Rcond+Rthppbulk+Rthppcnd;
%Rjunction = 0.03
Rtt1 = Rconv+Rcond2+Rthppcnd;%+Rjunction;
Rtt2 = Rconv+Rcond2+Rthppbulk+Rthppcnd;%+Rjunction;
rttt = Rtt1;
RSA = Rtt1+Rtt2;
dt1 = q_s*Rtt1;
dt2 = q_s*Rtt2+q_s/(massf*cp_liquid);
dtavg = (q_s*Rtt1+q_s*Rtt2+q_s/(massf*cp_liquid))/2;
dtf =q_s/(massf*cp_liquid);

dtfinal = dtavg+Tin;
massreq = massf*nm;
end