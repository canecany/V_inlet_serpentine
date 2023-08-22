function h = heat_transfer_coefficent(Pr,Re,k_f,Dh)
%function to find heat transfer coefficient 
% k_f = fluid thermal conductivity (W/m*K)
% Dh = hydralic diameter
if Re <= 2300 % Laminar flow
    Nu = 3.66;
elseif Re <= 4000 % Transition flow
    Nu = 0.023*Re^(4/5)*Pr^(1/3);
else % Turbulent flow
    Nu = 0.027*Re^(4/5)*Pr^(1/3);
end
h = Nu*k_f/Dh;
end