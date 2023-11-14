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
c       = 1;  % Airfoil chord
Uinf    = 1;   % Freestream Velocity field module
Naux       = [16,32,64,128,256,512];  % Number of panels
AoAaux     = [0,2,4,6,8,10];  % Angle of attack in degrees

% Vector definition
CL = zeros(size(Naux,2),size(AoAaux,2));
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
figure
hold on
for j=1:numel(AoAaux)
    AoA = AoAaux(j);
    Qinf = Uinf*[cosd(AoA);sind(AoA)];

    % POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
    % [gamma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N);
    [gamma,uInd,wInd] = computeConstantVortexDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N);
    
    % Preprocessing computations
    V   = computeVelocity(Qinf,gamma,uInd,wInd,N);
    cp  = computeCp(Qinf,V);
    CP(:,j) = cp;
    cl = computeCl(cp,lp,Ncj,c,AoA);
    CL(i,j) = cl;
    cm4 = computeCm4(cp,coord_xC,coord_xP,c);
    CM4(i,j) = cm4; 

    isentropicExp = 1.4;
    % Laitone's Rule
    Cp_0 = min(cp);
    syms Minf2
    eqn_Cp = Cp_0/(sqrt(1-Minf2)+(Cp_0*Minf2)/(2*sqrt(1-Minf2))*(1+Minf2*(isentropicExp-1)/2)) == 2/(isentropicExp*Minf2)*(((2+(isentropicExp-1)*Minf2)/(1+isentropicExp))^(isentropicExp/(isentropicExp-1))-1);
    Mcrit2 = solve(eqn_Cp, Minf2);
    Mcrit = double(sqrt(Mcrit2(1)));
    MCR(i,j) = Mcrit;

    msg = sprintf('Cl=%i and Cm1/4=%i for AoA=%i degrees', cl, cm4, AoA);
    disp(msg);

plot(coord_xC(:,1),cp);
end
title("Pressure coefficient plot for " + string(N) + " panels - $C_p~vs~x/c$")
xlabel("Chord length $x/c$ ");
ylabel("Pressure Coefficient $C_{p}$");
legend("$\alpha=0$","$\alpha=2$","$\alpha=4$","$\alpha=6$","$\alpha=8$","$\alpha=10$","Location","southeast");
grid on;
grid minor;
box on;
axis padded
fontsize(13,"points")
hold off;


% POSTPROCESSING
% plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj); % Panel and norm vector visualization 
% plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N);
% plotVelocityDistribution(Qinf,V,N);
% plotPressureCoefficient(coord_xP,coord_xC,Ncj,CP(:,a),N)
end

% MESH INDEPENDENCE TEST - PLOTS
figure
hold on
title("Panel density independence test - $C_l~vs~\alpha$")
plot(AoAaux,CL(1,:));
plot(AoAaux,CL(2,:));
plot(AoAaux,CL(3,:));
plot(AoAaux,CL(4,:));
plot(AoAaux,CL(5,:));
plot(AoAaux,CL(6,:));
xlabel("Angle of attack $\alpha$ ($^\circ$)");
ylabel("Lift Coefficient $C_{l}$");
legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","northwest");
grid on;
grid minor;
box on;
axis padded
fontsize(13,"points")
hold off;

figure
hold on
title("Panel density independence test - ${C_m}_{1/4}~vs~\alpha$")
plot(AoAaux,CM4(1,:));
plot(AoAaux,CM4(2,:));
plot(AoAaux,CM4(3,:));
plot(AoAaux,CM4(4,:));
plot(AoAaux,CM4(5,:));
plot(AoAaux,CM4(6,:));
xlabel("Angle of attack $\alpha$ ($^\circ$)");
ylabel("Momentum Coefficient about the quarter ${C_{m}}_{1/4}$");
legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","southwest");
grid on;
grid minor;
box on;
axis padded
fontsize(13,"points")
hold off;

figure
hold on
title("Panel density independence test - $M_{crit}~vs~\alpha$")
plot(AoAaux,MCR(1,:));
plot(AoAaux,MCR(2,:));
plot(AoAaux,MCR(3,:));
plot(AoAaux,MCR(4,:));
plot(AoAaux,MCR(5,:));
plot(AoAaux,MCR(6,:));
xlabel("Angle of attack $\alpha$ ($^\circ$)");
ylabel("Critical Freestream Mach number $M_{crit}$");
legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","southwest");
grid on;
grid minor;
box on;
axis padded
fontsize(13,"points")
hold off;


