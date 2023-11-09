%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; %close all;

% Input parameters
N             = 200;  % Number of panels
R             = 1;    % Radius of the cilinder
c             = 2*R;  % Equivalent chord
AoA           = 6;    % Angle of attack
isentropicExp = 1.4; % Gamma ideal gas
Mcrit         = 0.336947168481477;
a             = 340;  % [m/s] Sound velocity
M             = linspace(0.000001,Mcrit,25);
clCorrected   = zeros(size(M,2),1);
for i=1:size(M,2)
Uinf    = M(1,i)*a;   % Freestream Velocity field module
B       = sqrt(1-M(1,i)^2); %Prandlt-Glauert Correction Parameter
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
cpCorrected = cp./(sqrt(1-M(1,i)^2)+(cp.*M(1,i)^2)/(2*sqrt(1-M(1,i)^2))*(1+M(1,i)^2*(isentropicExp-1)/2));
cl = computeCl(cpCorrected,lp,Ncj,c,AoA);
clCorrected(i) = cl;

% POSTPROCESSING
% plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
% plotVelocityDistribution(Qinf,V,N);
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)

end

figure
plot(M,clCorrected);
