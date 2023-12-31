%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; %close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Input parameters
N             = 512;  % Number of panels
NACA          = "0010";
% R             = 1;    % Radius of the cilinder
c             = 1;  % Equivalent chord
AoA           = 2;    % Angle of attack
isentropicExp = 1.4; % Gamma ideal gas
Mcrit         = 0.604052986849087;
a             = 340;  % [m/s] Sound velocity
M             = (Mcrit-0.15):0.01:Mcrit;

% Vector definition
clCorrected   = zeros(size(M,2),1);

for i=1:size(M,2)
Uinf    = M(1,i)*a;   % Freestream Velocity field module
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
cp = computeCp(Qinf,V,gamma);
cpCorrected = cp./(sqrt(1-M(1,i)^2)+(cp.*M(1,i)^2)/(2*sqrt(1-M(1,i)^2))*(1+M(1,i)^2*(isentropicExp-1)/2));
[cl_int,cl_kutta] = computeCl(cpCorrected,lp,Ncj,c,AoA,Qinf,gamma);
clCorrected(i) = cl_int;

end

load RequestedDataFor3.mat;
figure
hold on
title("$C_l$ vs $M_\infty$")
plot(M,clCorrected,'b');
scatter(RequestedDataFor3.M,RequestedDataFor3.clCorrected,"square","filled","r");
xlabel("Freestream Mach number $M_{\infty}$");
ylabel("Lift Coefficient $C_{l}$");
grid on;
grid minor;
box on;
axis padded
hold off;

cpast = -0.7708; % Imposed Minimum Cp* to find Critical Mach number
