%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;

%% EXERCISE I
% Input parameters
N       = 200;  % Number of panels
R       = 1;    % Radius of the cilinder
AoA     = 6;    % Angle of attack
Uinf    = 1;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(AoA);sind(AoA)]; % Freestream Velocity field

% Precomputations
[coord_xP]      = setCylinderNodes(R,N); % Normally the coordinates xP are given
[coord_xC,lp]   = setGeometricParameters(coord_xP,N);
[cj,sj,Ncj,Tcj] = computePanelAngleAndNormalAndTangentVectors(coord_xP,lp,N); % Panel angle, normal and tangent vectors calculation

% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
% [gamma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N);
[gamma,uInd,wInd] = computeConstantVortexDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N);

% Preprocessing computations
V   = computeVelocity(Qinf,gamma,uInd,wInd,N);
cp  = computeCp(Qinf,V);

% POSTPROCESSING
plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
plotVelocityDistribution(Qinf,V,N);
plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)

%% EXERCISE II
isentropicExp = 1.4;
% Entropia
Cp_0 = min(cp);
syms Minf2
eqn_Cp = Cp_0/(sqrt(1-Minf2)+(Cp_0*Minf2)/(2*sqrt(1-Minf2))*(1+Minf2*(isentropicExp-1)/2)) == 2/(isentropicExp*Minf2)*(((2+(isentropicExp-1)*Minf2)/(1+isentropicExp))^(isentropicExp/(isentropicExp-1))-1);
Mcrit2 = solve(eqn_Cp, Minf2);
Mcrit = double(sqrt(Mcrit2(1)));


