function [rho, mu, nu, cp_liquid, lambda,Pr,L] = fluid_properties_code(Tin,P0,ele)
%function to fine the property of fluid
rho = py.CoolProp.CoolProp.PropsSI("D", "T", Tin, "P", P0, ele);
mu = py.CoolProp.CoolProp.PropsSI("V", "T",Tin , "P", P0, ele);
nu =mu/rho;
cp_liquid =py.CoolProp.CoolProp.PropsSI("C", "T", Tin , "P", P0, ele); %calculate Cp
Pr = py.CoolProp.CoolProp.PropsSI("PRANDTL", "T", Tin, "P", P0, ele);
L = py.CoolProp.CoolProp.PropsSI("L", "T", Tin, "P", P0, ele);
lambda = nu/(Pr/(rho*cp_liquid));
end