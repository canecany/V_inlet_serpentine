function [dtfinal,vreq,massreq,nm,Tr,RSA,Re] = heatf(Tmax,Tin,W,L,tc,qc,td,nt)
%input variable
% Tmax = average temperature of The cold plate design  degC
% Tin = temperature inlet degC
% W = width of PEEBS m
% L = Length of pebbs m
% tc = thickness oc cool plate %inch
% qc = heat form Pebbs W
% td = tube diameter %inch
% nt = number of turn inthe cold plate

%output variable
% dt final = temperature average of plate
% vreq = velocity inlet require use for determine the function of the presure drop
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
vmin = Low_Re*mu/(rho*Dm);
n =1000;
vin = linspace(vmin,2.5,n);      %velocity inlet per chanel
for i =1 :n
    Re = rho*vin(i)*Dm/mu                                             %correct
    rei(i) = Re;                                                        %correct
    h = heat_transfer_coefficent(Pr,Re,Cond_f,Dm);                     %correct
    massf = rho*Acss*vin(i);
    mass_f_plot(i) =massf;
    ht(i) = h;
    Rconv = 1/(h*HTA);
    Rthppbulk = (L-Tr*(nt/2))*(nt+1)*(2*pi*Dm)/(massf*cp_liquid);     %neet to check this a bit
    Rcond1 = ((pi*Dm)/(pi*k_s)) * log (1/(sin(0.25*pi*Dm/(pi*Dm))));
    S = 2*pi*(L-Tr*(nt/2))*(nt+1)/(log(1.08*t/Dm));
    Rcond2 = 1/(S*k_s);
    Rthppcnd = tsub/k_s;
    %Rtt = Rconv+Rcond+Rthppbulk+Rthppcnd;
    %Rjunction = 0.03
    Rtt1 = Rconv+Rcond2+Rthppcnd;%+Rjunction;
    Rtt2 = Rconv+Rcond2+Rthppbulk+Rthppcnd;%+Rjunction;
    rttt(i) = Rtt1;
    RSA = Rtt1+Rtt2;
    dt1(i) = q_s*Rtt1;
    dt2(i) = q_s*Rtt2+q_s/(massf*cp_liquid);
    dtavg(i) = (q_s*Rtt1+q_s*Rtt2+q_s/(massf*cp_liquid))/2;
    dtf(i) =q_s/(massf*cp_liquid);
    if dtavg(i) < Tmax-Tin
        dtfinal = dtavg(i)+Tin;
        massreq = massf*nm;
        vreq = vin(i);
        disp(dtfinal)
        break
    end

    if i == n
        dtfinal = 0;
        massreq = 0;
        vreq = 0;
        nm = 0;
        disp('unrealistic value')
    end
    
end
end