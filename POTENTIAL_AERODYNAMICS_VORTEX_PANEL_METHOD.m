%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% EXERCISE I
% Input parameters
N       = 512;  % Number of panels
NACA    = "0010";
% R       = 1;    % Radius of the cilinder
c       = 1;  % Airfoil chord
AoA     = 2;  % Angle of attack
Uinf    = 1;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(AoA);sind(AoA)]; % Freestream Velocity field

% Precomputations

% Cilinder
% [coord_xP]      = setCylinderNodes(R,N); % Normally the coordinates xP are given
% [coord_xC,lp]   = setGeometricParameters(coord_xP,N);

% Airfoil
[coord_xP,coord_xC,lp] = setGeometricParameters(c,N,NACA);
[cj,sj,Ncj,Tcj] = computePanelAngleAndNormalAndTangentVectors(coord_xP,lp,N); % Panel angle, normal and tangent vectors calculation

% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
% [gamma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N);
[gamma,uInd,wInd] = computeConstantVortexDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N);

% Preprocessing computations
V   = computeVelocity(Qinf,gamma,uInd,wInd,N);
cp  = computeCp(Qinf,V);
cl = computeCl(cp,lp,Ncj,c,AoA);
cm4 = computeCm4(cp,coord_xC,coord_xP,c);
% POSTPROCESSING
plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
plotVelocityDistribution(Qinf,V,N);
plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)

%% EXERCISE II
isentropicExp = 1.4;
% Laitone's Rule
Cp_0 = min(cp);
syms Minf2
eqn_Cp = Cp_0/(sqrt(1-Minf2)+(Cp_0*Minf2)/(2*sqrt(1-Minf2))*(1+Minf2*(isentropicExp-1)/2)) == 2/(isentropicExp*Minf2)*(((2+(isentropicExp-1)*Minf2)/(1+isentropicExp))^(isentropicExp/(isentropicExp-1))-1);
Mcrit2 = solve(eqn_Cp, Minf2);
Mcrit = double(sqrt(Mcrit2(1)));


