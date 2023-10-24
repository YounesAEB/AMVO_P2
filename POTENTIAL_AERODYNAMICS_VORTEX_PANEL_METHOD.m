%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;

% Input parameters
N       = 50;  % Number of panels
R       = 1;    % Radius of the cilinder
AoA     = 6;    % Angle of attack
Uinf    = 30;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(AoA);sind(AoA)]; % Freestream Velocity field

% Precomputations
[coord_xP]      = setCylinderNodes(R,N); % Normally the coordinates xP are given
[coord_xC,lp]   = setGeometricParameters(coord_xP,N);
[cj,sj,Ncj,Tcj] = computePanelAngleAndNormalAndTangentVectors(coord_xP,lp,N); % Panel angle, normal and tangent vectors calculation

% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
%[gamma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N);
[gamma,uInd,wInd] = computeConstantVortexDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N);

% Preprocessing computations
V   = computeVelocity(Qinf,gamma,uInd,wInd,N);
cp  = computeCp(Qinf,V);


% POSTPROCESSING
plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
plotVelocityDistribution(Qinf,V,N);
plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N);


