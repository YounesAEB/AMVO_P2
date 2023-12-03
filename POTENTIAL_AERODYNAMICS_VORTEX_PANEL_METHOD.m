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
NACA    = "0010";
c       = 1;   % Airfoil chord
Uinf    = 1;   % Freestream Velocity field module

% DO NOT FORGET TO CHANGE THE FOLLOWING INPUTS DEPENDING ON THE REQUIRED PLOT
Naux    = [512]; % Number of panels
AoAaux  = 0:0.5:6;  % Angle of attack in degrees

% Vector definition
CL_int = zeros(size(Naux,2),size(AoAaux,2));
CL_kutta = zeros(size(Naux,2),size(AoAaux,2));
CM4 = zeros(size(Naux,2),size(AoAaux,2));
MCR = zeros(size(Naux,2),size(AoAaux,2));

for i=1:size(Naux,2)
% R     = 1;    % Radius of the cilinder
N = Naux(i);

msg = sprintf('Calculation for %i panels...', N);
    disp(msg);

% Precomputations
% [coord_xP]      = setCylinderNodes(R,N); % Normally the coordinates xP are given
[coord_xP,coord_xC,lp] = setGeometricParameters(c,N,NACA);
[cj,sj,Ncj,Tcj] = computePanelAngleAndNormalAndTangentVectors(coord_xP,lp,N); % Panel angle, normal and tangent vectors calculation

% Vector definition
CP = zeros(N,size(AoAaux,2));

for j=1:numel(AoAaux)
    AoA = AoAaux(j);
    Qinf = Uinf*[cosd(AoA);sind(AoA)];

    % POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
    % [gamma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N);
    [gamma,uInd,wInd] = computeConstantVortexDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N);
    
    % Preprocessing computations
    V                   = computeVelocity(Qinf,gamma,uInd,wInd,N);

    cp                  = computeCp(Qinf,V,gamma); %Pressure Coefficient
    CP(:,j)             = cp;

    [cl_int,cl_kutta]   = computeCl(cp,lp,Ncj,c,AoA,Qinf,gamma); %Lift coefficient
    CL_int(i,j)         = cl_int;
    CL_kutta(i,j)       = cl_kutta;

    cm4                 = computeCm4(cp,coord_xC,coord_xP,c); %1/4 pitching moment coefficient
    CM4(i,j)            = cm4; 

    Mcrit = computeCriticalMach(cp);
    MCR(i,j) = Mcrit;

    msg = sprintf('Cl=%i for the integral, Cl=%i for kutta and Cm1/4=%i for AoA=%i degrees', cl_int, cl_kutta, cm4, AoA);
    disp(msg);

    %plotCoefficientPressure(coord_xC,cp,N,j,AoAaux);
end


% POSTPROCESSING
% plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
% plotVelocityDistribution(Qinf,V,N);
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,CP(:,a),N)
end

% PLOTS
%plotMeshIndepenceTests(AoAaux,Naux,CL_int,CM4,MCR);
plotExperimentalAndPanelMethodCl(AoAaux,CL_kutta);

% CODE TO PRINT THE FIGURES IN PDF FORMAT
%     set(gcf, 'Units', 'Centimeters');
%     pos = get(gcf, 'Position');
%     set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%         'PaperSize',[pos(3), pos(4)]);
%     print(gcf, 'Name', '-dpdf', '-r0'); % incrementar '-r0' resolución 

%% POST-PROCESS FOR PART 2

x  = AoAaux(:); %reshape the data into a column vector
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
plot(x0,y0,'gx','linewidth',4)
% Plot fitted data
plot(x,yhat,'r','linewidth',2)
hold off
