% MATLAB R2019b
% Dated: July-28-2020
% Author: Chi Chen @ Boston University
% Email: chenchi@bu.edu

% This script provides the partial derivatives of Ts with respect to albedo, 
% aerodynamic resistance and surface resistance. The Ts is analytically 
% solved with first-order Taylor series expansion and second-order Taylor 
% series expansion based on the land surface energy balance equation.

% Academic use ONLY, and please acknowledge the following paper:
% Blank for the publication
% Blank for the publication

clear
clc

%% Definitions
% FOTSE    first-order Taylor series expansion
% SOTSE    second-order Taylor series expansion
% Ts_1st   Ts solved with first-order Taylor series expansion
% Ts_2nd   Ts solved with second-order Taylor series expansion

%% constants
% Cp      specific heat of air at constant pressure  1004.64 J kg^-1 K^(-1)
% Lv      latent heat vaporization                    2.4665x10^6 J kg^(-1)
%      or latent heat of sublimation                  2.8002x10^6 J kg^(-1)
% sigma   Stefan-Boltzmann constant        5.670367x10^(-8) W m^(-2) K^(-4)
% Rwa     Ratio molecular weight of water vapor/dry air               0.622

%% parameters
% albedo  albedo                                              dimensionless
% emis    surface emissivity                                  dimensionless
% e_sat   saturation vapor pressure                                      Pa
% G       ground heat flux                                         W m^(-2)
% Lin     incoming longwave radiation                              W m^(-2)
% P       atmospheric pressure                                           Pa
% qa      air specific humidity                                   g kg-^(1)
% q_sat   saturation specific humidity                            g kg-^(1)
% rhoa    air density                                             kg m^(-3)
% ra      aerodynamic resistance                                   s m^(-1)
% rs      surface resistance                                       s m^(-1)
% Sin     incoming shortwave radiation                             W m^(-2)
% Ta      air temperature above blending height                           K
% Ts      land surface temperature                                        K

%% intemediate parameters / shortcuts
% a
% b
% c
% delta
% gamma
% f
% lambda0prime
% r0
% Rn_star

syms Cp Lv sigma Rwa% constants
syms albedo emis e_sat G Lin P qa q_sat rhoa ra rs Sin Ta Ts % parameters
syms delta gamma f lambda0prime r0 Rn_star a b c % shortcuts

%% FOTSE LST model, i.e., Ts_1st
% In the paper: e_sat = 611*exp(17.27*(273.15-Ta)/(Ta-35.85)); % (Dingman, 2008)
q_sat = Rwa / P * e_sat;

% define shortcuts
delta        = diff(e_sat, Ta,1);
gamma        = (Cp*P)/(Rwa*Lv);
lambda0prime = 1/(rhoa*Cp);
r0           = rhoa*Cp/(4*emis*sigma*Ta^3);
Rn_star      = Sin*(1-albedo) + emis*Lin -emis*sigma*Ta^4;
f            = 1/r0+(1/ra)*(1 + (delta/gamma)*(ra/(ra + rs)));

% analytical expression for Ts_1st
Ts_1st       = lambda0prime*(Rn_star-G-rhoa*Lv*(q_sat-qa)/(ra+rs))/f + Ta;

% 1st order derivatives for FOTSE LST model
dTs_1st_dalbedo  = diff(Ts_1st, albedo, 1);
dTs_1st_dra      = diff(Ts_1st, ra, 1);
dTs_1st_drs      = diff(Ts_1st, rs, 1);

% 2nd order derivatives for FOTSE LST model
ddTs_1st_ddalbedo = diff(Ts_1st, albedo, 2);
ddTs_1st_ddra     = diff(Ts_1st, ra, 2);
ddTs_1st_ddrs     = diff(Ts_1st, rs, 2);

% degree 2 cross-order derivatives for FOTSE LST model
ddTs_1st_dalbedo_dra = diff(diff(Ts_1st, albedo, 1), ra, 1);
ddTs_1st_dra_drs     = diff(diff(Ts_1st, ra, 2), rs, 1);
ddTs_1st_drs_dalebdo = diff(diff(Ts_1st, rs, 1), albedo, 1);

%% SOTSE LST model, i.e., Ts_2nd
% shortcuts
a = 6*emis*sigma*Ta^2 + 1/2 * diff(q_sat,Ta,2);
b = 4*emis*sigma*Ta^3 + rhoa*Lv*diff(q_sat,Ta,1)/(ra+rs) + rhoa*Cp/ra;
c = rhoa*Lv*(q_sat-qa)/(ra+rs) - (Rn_star-G);

% analytical expression for Ts_2nd
Ts_2nd = -b + sqrt(b^2-4*a*c) / (2*a) + Ta;

% 1st order derivatives for SOTSE LST model
dTs_2nd_dalbedo  = diff(Ts_2nd, albedo, 1);
dTs_2nd_dra      = diff(Ts_2nd, ra, 1);
dTs_2nd_drs      = diff(Ts_2nd, rs, 1);

% 2nd order derivatives for FOTSE LST model
ddTs_2nd_ddalbedo = diff(Ts_2nd, albedo, 2);
ddTs_2nd_ddra     = diff(Ts_2nd, ra, 2);
ddTs_2nd_ddrs     = diff(Ts_2nd, rs, 2);

% degree 2 cross-order derivatives for FOTSE LST model
ddTs_2nd_dalbedo_dra = diff(diff(Ts_2nd, albedo, 1), ra, 1);
ddTs_2nd_dra_drs     = diff(diff(Ts_2nd, ra, 2), rs, 1);
ddTs_2nd_drs_dalebdo = diff(diff(Ts_2nd, rs, 1), albedo, 1);

% END
