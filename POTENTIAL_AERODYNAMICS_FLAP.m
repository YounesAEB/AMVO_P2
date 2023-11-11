%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;

%% COMPOUND WING
% Input parameters
N       = 512;  % Number of panels main airfoil
M       = 512;   % Number of panels flap airfoil
NACA    = "0015";
c1      = 1;  % Main airfoil chord
c2      = 0.45;  % Flap airfoil chord
d       = 0.05; % Gap
c       = c1 + c2 + d; % Total chord
AoA     = 0;  % Angle of attack main airfoil
Uinf    = 1;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(AoA);sind(AoA)]; % Freestream Velocity field

df_aux      = [4,8,12,16,20];  % Flap deflection
CM4 = zeros(size(df_aux,2),1);
CL  = zeros(size(df_aux,2),1);
CL1  = zeros(size(df_aux,2),1);
CL2  = zeros(size(df_aux,2),1);
for i=1:size(df_aux,2)

df = df_aux(i);

% Precomputations
[coord_xP,coord_xC,lp] = setFlapGeometry(c1,c2,d,df,N,M,NACA);
[cjN,sjN,NcjN,TcjN] = computePanelAngleAndNormalAndTangentVectors(coord_xP(1:N+1,:),lp(1:N,:),N); % Main panel angle, normal and tangent vectors calculation
[cjM,sjM,NcjM,TcjM] = computePanelAngleAndNormalAndTangentVectors(coord_xP(N+2:end,:),lp(N+1:end,:),M); % Flap panel angle, normal and tangent vectors calculation
cj = [cjN;cjM];
sj = [sjN;sjM];
Ncj = [NcjN;NcjM];
Tcj = [TcjN;TcjM];

% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
[gamma,uInd,wInd] = computeConstantVortexDistributionFlap(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N,M);

% Preprocessing computations
V   = computeVelocity(Qinf,gamma,uInd,wInd,N+M);
cp  = computeCp(Qinf,V);
cl1 = computeCl(cp(1:N,:),lp(1:N,:),Ncj(1:N,:),c,AoA);
cl2 = computeCl(cp(N+1:N+M,:),lp(N+1:N+M,:),Ncj(N+1:N+M,:),c,AoA);
cl = cl1 + cl2;
cm4 = computeCm4Flap(cp,coord_xC,coord_xP,c,N,M);
CM4(i,1) = cm4;
CL(i,1) = cl;
CL1(i,1) = cl1;
CL2(i,1) = cl2;
end

% POSTPROCESSING
% plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
% plotVelocityDistribution(Qinf,V,N);
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)


