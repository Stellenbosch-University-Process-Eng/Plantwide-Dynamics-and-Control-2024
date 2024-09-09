function [FST, dBx_dt, dTJ_dt, dP_dt, dLJ_dt] = MultiEffectEvaporator_ODEs(x, P_upstream,Vdraw_total, FJ_in, Bx_in, TJ_in, FJ_out, Bx, TJ, P, LJ)
% The model is based on the MEng theses by PD Smith (2000)
% "Control and optimization of a multiple-effect evaporator"
% University of Capte Town, http://hdl.handle.net/11427/5397
%
% Inputs:
%  FJ_in (t/hr): mass flowrate of juice entering the evaporator
%  FJ_out (t/hr): mass flowrate of juice leaving the evaporator
%  Bx_in (%): brix of juice entering the evaporator
%  TJ_in (oC): temperature of juice entering the evaporator
%  TST (oC): temperature of incoming steam
%  Vnext (t/hr): mass flowrate of vapour to next effect
%  Vdraw (t/hr): mass flowrate of vapour drawn to equipment elsewhere
%
% State variables:
%  P (kPa): vapour phase absolute pressure
%  mJ (t): mass of juice
%  Bx (%): brix of juice
%  TJ (oC): temperature of juice

%% Parameters
R = 8.314*(1e3/18);   % (m3.kPa/t.K): gas constant for water
Vol = 120;            % (m3): vapour space volume
alpha = 1e2;          % (t/hr.kPa): coefficient relating vapourisation rate to pressure driving force
U = 1500*3600*1e-6;   % (W/K.m2)*3600(s/hr)*1e-6(MJ/J) = (MJ/hr.K.m2): base heat transfer coefficient

Dves = 5.5;           % (m): vessel diameter
CrossA = pi*Dves^2/4;    % (m2): juice cross-sectional area

Ntubes = 7550;        % (~): number of tubes for heat transfer
Ltubes = 2;           % (m): length of tube
Dtubes = 42e-3;       % (m): tube diameter
A = Ntubes*pi*Dtubes*Ltubes; % (m2): heat transfer area 

%% Mass of juice
mJ = LJ*frho(TJ, Bx)*CrossA; % (t): mass of juice

%% Vapour evolution
% Driving force for vapour evolution: sugar decreases the vapour pressure
% exerted by the juice, hence the driving force for vapourisation
Vev = alpha*(f_PJ_eq(TJ, Bx) - P); % (t/hr): mass flowrate of vapour evolved from boiling juice

%% Heat transfer 
% We assume the steam is saturated
PST = (x)*P_upstream + (1-x)*fP_eq(TJ);
TST = fT_eq(PST);
Q = U*A*(TST - TJ);          % (MJ/hr): heat transferred from steam

% Steam flowrate required to maintain heat transfer rate
FST = Q/fdhvap(TST);              % (t/hr): steam flowrate

%% Differential equations
% Mass balance, and rate of change of level
dmJ_dt = ( FJ_in - FJ_out - Vev );  
dLJ_dt = dmJ_dt/(frho(TJ, Bx)*CrossA);       

% Solids balance
dBx_dt = ( FJ_in*Bx_in - FJ_out*Bx - Bx*dmJ_dt ) / mJ;

% Energy balance
dTJ_dt =  ( FJ_in*fh(TJ_in, Bx_in) - FJ_out*fh(TJ, Bx) - Vev*fh_ev(TJ, Bx) + Q ) / ( mJ*fcP(TJ,Bx) )...
        - TJ*dmJ_dt/mJ; 

% Vapour phase mass balance
dP_dt =  ( R*TJ/Vol )*( Vev - Vdraw_total ) + P/TJ*dTJ_dt;

end

%% Thermodynamic relationships
% Saturation pressure (pure water)
function P_eq = fP_eq(T)
    P_eq = 10^( 7.8656 - 2188.8/(T+273.15) );
end

% Saturation temperature (pure water)
function T_eq = fT_eq(P)
    T_eq =2188.8/(7.8656 - log10(P)) - 273.15;
end

% Saturation pressure (juice)
function P_eq = f_PJ_eq(T, Bx)
    BPE = 6.054e-5*(  ( (T+273.15)^2*Bx^2 )/( (374.3-T)^2 ) * 5.84e-7*(Bx-40)^2 + 7.2e-4); % (oC):  boiling point elevation
    Tvap = T - BPE;        % (oC):  equilibrium temperature of vapour evolved from boiling juice
    P_eq = fP_eq(Tvap); % (kPa): equilibrium pressure of vapour evolved from boiling juice
end

% Specific enthalpy (juice)
function h = fh(T,Bx)
    h = 2.326*( (Bx/10)*(100+Bx)/(900 - 8*Bx) ...
                         + 1.8*T*(1-Bx/100)*(0.6-9e-4*T) ); % (MJ/t), enthalpy of juice
end

% Heat of vaporisation (pure water)
function dhvap = fdhvap(T)
    Tr = (T + 273.15)/647.13;  % (~): reduced temperature of outgoing steam
    dhvap = 2889*(1-Tr)^(0.3199 - 0.212*Tr + 0.25795*Tr^2); % (MJ/t): heat of vapourisation
end

% Specific enthalpy (water vapour)
function hev = fh_ev(T,Bx)
    % Assuming thermal equilibrium, the temperature of the steam going out is
    % equal to the temperature of the juice (slightly overheated). This informs the specific heat
    hev = fh(T, Bx) + fdhvap(T); % (MJ/t): specific enthalpy of vapour evolved from boiling juice
end

% Specific heat capacity (juice)
function cP = fcP(T, Bx)
    cP =   4.1253 - 0.02804*Bx + 6.7e-5*Bx*T...
         + 1.8691e-3*T - 9.271e-6*T^2; % (MJ/t.K): specific heat capacity of the juice
end

% Density (juice)
function rho = frho(T, Bx)
    rho = ( 1 + Bx*(Bx + 200)/54e3 ) * ( 1 - 0.036*(T-20)/(160-T) );     % (t/m3): juice density
end