%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;

%% COMPOUND WING
% Input parameters
N       = 32;  % Number of panels main airfoil
M       = 32;   % Number of panels flap airfoil
NACA    = "0015";
c1      = 1;  % Main airfoil chord
c2      = 0.45;  % Flap airfoil chord
d       = 0.05; % Gap
c       = c1 + c2 + d; % Total chord
AoA     = 0;  % Angle of attack main airfoil
df      = 4;  % Flap deflection
Uinf    = 1;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(AoA);sind(AoA)]; % Freestream Velocity field

% Precomputations
[coord_xP,coord_xC,lp] = setFlapGeometry(c1,c2,d,df,N,M,NACA);
[cjN,sjN,NcjN,TcjN] = computePanelAngleAndNormalAndTangentVectors(coord_xP(1:N+1,:),lp(1:N,:),N); % Main panel angle, normal and tangent vectors calculation
[cjM,sjM,NcjM,TcjM] = computePanelAngleAndNormalAndTangentVectors(coord_xP(N+2:N+M+2,:),lp(N+1:N+M,:),M); % Flap panel angle, normal and tangent vectors calculation
cj = [cjN;cjM];
sj = [sjN;sjM];
Ncj = [NcjN;NcjM];
Tcj = [TcjN;TcjM];

% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
[gamma,uInd,wInd] = computeConstantVortexDistributionFlap(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N,M);

% % Preprocessing computations
% V   = computeVelocity(Qinf,gamma,uInd,wInd,N);
% cp  = computeCp(Qinf,V);
% cm4 = computeCm4(cp,coord_xC,coord_xP,c);
% POSTPROCESSING
plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
% plotVelocityDistribution(Qinf,V,N);
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)


