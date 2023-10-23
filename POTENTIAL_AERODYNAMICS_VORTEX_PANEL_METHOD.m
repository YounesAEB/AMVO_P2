%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - AMVO 
%  /  Matlab code to assess the numerical solution of potential equations                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;

% Input parameters

N = 80; % Number of panels
R = 1; % Radius of the cilinder
Qinf = [30;0]; % Freestream Velocity field

%% GEOMETRIC PREPROCESS - THIS CAN ALL BE A FUNCTION!
% Coordinate definition and calculation
% Solved for a 2D cilinder, for complex geometries compute xC, xP are given

coord_xP = zeros(N+1,2); % Panel limits coordinates
coord_xC = zeros(N,2); % Panel center nodes

for i = 1:N+1
    coord_xP(i,1) = R*cosd(-(360/N)*(i-1));
    coord_xP(i,2) = R*sind(-(360/N)*(i-1));
end

for i = 1:N
    coord_xC(i,1) = (coord_xP(i,1) + coord_xP(i+1,1))/2;
    coord_xC(i,2) = (coord_xP(i,2) + coord_xP(i+1,2))/2;
end

% Panel length definition and calculation

lp = zeros(N,1);

for i = 1:N
    lp(i) = sqrt((coord_xP(i+1,1)-coord_xP(i,1))^2 + ...
                 (coord_xP(i+1,2)-coord_xP(i,2))^2);
end

% Panel angle calculation and normal and tangent vectors calculation

cj  = zeros(N,1);
sj  = zeros(N,1);
Ncj = zeros(N,2);
Tcj = zeros(N,2);

for i = 1:N
    cj(i) = (coord_xP(i+1,1)-coord_xP(i,1))/lp(i);
    sj(i) = (coord_xP(i,2)-coord_xP(i+1,2))/lp(i);
    Ncj(i,:)   = [sj(i) , cj(i)];
    Tcj(i,:)   = [cj(i) , -sj(i)];
end

% Panel and norm vector visualization - This can be an independent function
figure()
scatter(coord_xP(:,1),coord_xP(:,2),30,'k','filled');
hold on
plot(coord_xP(:,1),coord_xP(:,2),'k');
quiver(coord_xC(:,1),coord_xC(:,2),Ncj(:,1),Ncj(:,2),'r');
scatter(coord_xC(:,1),coord_xC(:,2),5,'r');
hold off
axis equal

%% POTENTIAL AERODYNAMICS - VELOCITY AND PRESSURE FIELDS CALCULATION
%%{
% CONSTANT STRENGHT SOURCE DISTRIBUTION - THIS CAN ALL BE A FUNCTION! 
% (Used for symmetric airfoils with zero angle of attack)
% Variable Definition
b = zeros(N,1); % Vector of independent terms
A = zeros(N,N); % Influence coefficient matrix 
xLoc = zeros(N,N); % X coordinate of the control point i in the local frame j       
zLoc = zeros(N,N); % Y coordinate of the control point i in the local frame j
uIndLoc = zeros(N,N); % Induced horizontal velocity at point i by panel j in the local frame j 
wIndLoc = zeros(N,N); % Induced vertical velocity at point i by panel j in the local frame j 

uInd = zeros(N,N); % Induced horizontal velocity at point i by panel j in the global frame
wInd = zeros(N,N); % Induced vertical velocity at point i by panel j in the global frame

for i = 1:N
    b(i) = -Ncj(i,:)*Qinf;
    for j = 1:N
        if i == j
            A(i,i) = 0.5;
        else
            xLoc(i,j) = (coord_xC(i,1)-coord_xP(j,1))*cj(j)...
                      - (coord_xC(i,2)-coord_xP(j,2))*sj(j);

            zLoc(i,j) = (coord_xC(i,1)-coord_xP(j,1))*sj(j)...
                      + (coord_xC(i,2)-coord_xP(j,2))*cj(j);
            
            r1 = sqrt((coord_xP(j,1)-coord_xC(i,1))^2 + ...
                       (coord_xP(j,2)-coord_xC(i,2))^2);

            r2 = sqrt((coord_xP(j+1,1)-coord_xC(i,1))^2 + ...
                       (coord_xP(j+1,2)-coord_xC(i,2))^2); 

            theta1 = atan2(zLoc(i,j),xLoc(i,j));

            theta2 = atan2(zLoc(i,j),xLoc(i,j)-lp(j));

            uIndLoc(i,j) = (1/4*pi)*log((r1^2)/(r2^2));
            wIndLoc(i,j) = (theta2 - theta1)/(2*pi);
            
            % Change of the velocity components to the global frame
            uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
            wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);

            A(i,j) = [uInd(i,j) wInd(i,j)]*Ncj(i,:)';
        end        
    end
end

% Source strenghts obtention 
gamma = A\b; % sigma in reality, stated gamma to use the same post process code
%}
%{
% CONSTANT STRENGHT VORTEX DISTRIBUTION - THIS CAN ALL BE A FUNCTION! 
% (Used for asymmetric or symmetric airfoils with non-zero angle of attack)
% Variable Definition
b = zeros(N,1); % Vector of independent terms
A = zeros(N,N); % Influence coefficient matrix 
xLoc = zeros(N,N); % X coordinate of the control point i in the local frame j       
zLoc = zeros(N,N); % Y coordinate of the control point i in the local frame j
uIndLoc = zeros(N,N); % Induced horizontal velocity at point i by panel j in the local frame j 
wIndLoc = zeros(N,N); % Induced vertical velocity at point i by panel j in the local frame j 

uInd = zeros(N,N); % Induced horizontal velocity at point i by panel j in the global frame
wInd = zeros(N,N); % Induced vertical velocity at point i by panel j in the global frame

for i = 1:N
    b(i) = -Tcj(i,:)*Qinf;
    for j = 1:N
        if i == j
            A(i,i) = -0.5;
        else
            xLoc(i,j) = (coord_xC(i,1)-coord_xP(j,1))*cj(j)...
                      - (coord_xC(i,2)-coord_xP(j,2))*sj(j);

            zLoc(i,j) = (coord_xC(i,1)-coord_xP(j,1))*sj(j)...
                      + (coord_xC(i,2)-coord_xP(j,2))*cj(j);
            
            r1 = sqrt((coord_xP(j,1)-coord_xC(i,1))^2 + ...
                       (coord_xP(j,2)-coord_xC(i,2))^2);

            r2 = sqrt((coord_xP(j+1,1)-coord_xC(i,1))^2 + ...
                       (coord_xP(j+1,2)-coord_xC(i,2))^2); 

            theta1 = atan2(zLoc(i,j),xLoc(i,j));

            theta2 = atan2(zLoc(i,j),xLoc(i,j)-lp(j));

            uIndLoc(i,j) = (theta2 - theta1)/(2*pi);
            wIndLoc(i,j) = (1/4*pi)*log((r2^2)/(r1^2));
            
            % Change of the velocity components to the global frame
            uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
            wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);

            A(i,j) = [uInd(i,j) wInd(i,j)]*Tcj(i,:)';
        end        
    end
end
% Kutta condition application (gamma_1 + gamma_N = 0)
A(N/4,:)   = 0;
A(N/4,1)   = 1;
A(N/4,end) = 1;
b(N/4)     = 0;
% Sorce strenghts obtention 
gamma = A\b;
%}
%% POSTPROCESSING
%%{
% Control point velocity computation: Vi = Qinf + sum_j=1^N(sig_j*V^~_i,j)
V = zeros(N,2);
for i=1:N
sum_i = 0;
sum_j = 0;
    for j=1:N
        sum_i = sum_i + gamma(j)*uInd(i,j);
        sum_j = sum_j + gamma(j)*wInd(i,j);
    end
V(i,1) = Qinf(1) + sum_i;
V(i,2) = Qinf(2) + sum_j;
end

% Control point pressure coefficient computation: Cp_i=1-((abs(Vi)/Qinf)^2
Cp = 1 - (sqrt(V(:,1).^2+V(:,2).^2)/sqrt(Qinf(1,1)^2+Qinf(2,1)^2)).^2;

% Source strenght distribution plot
intensity = zeros(N,2);
f         = 0.05; %Scale Factor
figure
for i = 1:N
    intensity(i,:) = f*gamma(i)*Ncj(i,:);
    if gamma(i) < 0
        quiver(coord_xC(i,1),coord_xC(i,2),-intensity(i,1),-intensity(i,2),'r');
        hold on
    else
        quiver(coord_xC(i,1),coord_xC(i,2),intensity(i,1),intensity(i,2),'b');
        hold on 
    end
end
plot(coord_xP(:,1),coord_xP(:,2),'k');
axis equal
hold off

% Pressure Coefficient chord distribution plot
figure
plot(coord_xC(:,1),Cp);

% Pressure Coefficient airfoil distribution plot
figure
for i = 1:N
    Cp_norm = Cp.*Ncj;
    if Cp(i) < 0
        quiver(coord_xC(i,1),coord_xC(i,2),-Cp_norm(i,1),-Cp_norm(i,2),'r');
        hold on
    else
        quiver(coord_xC(i,1),coord_xC(i,2),Cp_norm(i,1),Cp_norm(i,2),'b');
        hold on 
    end
end
plot(coord_xP(:,1),coord_xP(:,2),'k');
axis equal
hold off
% Streamlines arround airfoil plot

%}
