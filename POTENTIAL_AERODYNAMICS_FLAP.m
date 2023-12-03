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

df_aux      = 0:1:20;  % Flap deflection
% df_aux      = [4];  % Flap deflection
CM4 = zeros(size(df_aux,2),1);
CL_int  = zeros(size(df_aux,2),1);
CL1_int  = zeros(size(df_aux,2),1);
CL2_int  = zeros(size(df_aux,2),1);
CL_kutta  = zeros(size(df_aux,2),1);
CL1_kutta  = zeros(size(df_aux,2),1);
CL2_kutta  = zeros(size(df_aux,2),1);
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
cp = computeCp(Qinf,V,gamma);
[cl1_int,cl1_kutta] = computeCl(cp(1:N,:),lp(1:N,:),Ncj(1:N,:),c,AoA,Qinf,gamma(1:N,:));
[cl2_int,cl2_kutta] = computeCl(cp(N+1:N+M,:),lp(N+1:N+M,:),Ncj(N+1:N+M,:),c,AoA,Qinf,gamma(N+1:N+M,:));
cl_int = cl1_int + cl2_int;
cl_kutta = cl1_kutta + cl2_kutta;
cm4 = computeCm4Flap(cp,coord_xC,coord_xP,c,N,M);
CM4(i,1) = cm4;
CL_int(i,1) = cl_int;
CL1_int(i,1) = cl1_int;
CL2_int(i,1) = cl2_int;
CL_kutta(i,1) = cl_kutta;
CL1_kutta(i,1) = cl1_kutta;
CL2_kutta(i,1) = cl2_kutta;
    msg = sprintf('Cl=%i for the integral, Cl=%i for kutta and Cm1/4=%i for df=%i degrees', cl_int, cl_kutta, cm4, df);
    disp(msg);
end

% POSTPROCESSING
% plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N+M);
% plotVelocityDistribution(Qinf,V(1:512,:),N); % Main airfoil Velocity distribution
% plotVelocityDistribution(Qinf,V(512+1:end,:),M); % Flap airfoil Velocity distribution
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,M+N);

%% POST-PROCESS FOR PART 2

x  = df_aux(:); %reshape the data into a column vector
y  = CL_int(:);
x0 = 0;         %point to go through
y0 = 0;


% 'C' is the Vandermonde matrix for 'x'
n = 1; % Degree of polynomial to fit
Van(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    Van(:,j) = x.*Van(:,j+1);
end
C = Van;

% 'd' is the vector of target values, 'y'.
d = y;

% There are no inequality constraints in this case, i.e.,
A = [];
b = [];

% We use linear equality constraints to force the curve to hit the required point. In
% this case, 'Aeq' is the Vandermoonde matrix for 'x0'
Aeq = x0.^(n:-1:0);
% and 'beq' is the value the curve should take at that point
beq = y0;

p = lsqlin( C, d, A, b, Aeq, beq );

% We can then use POLYVAL to evaluate the fitted curve
yhat = polyval( p, x );

% Plot original data
plot(x,y,'.b-')
hold on
% Plot point to go through
plot(x0,y0,'gx','linewidth',2)
% Plot fitted data
plot(x,yhat,'r','linewidth',1)
hold off
