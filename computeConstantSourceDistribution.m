function [sigma,uInd,wInd] = computeConstantSourceDistribution(Qinf,coord_xP,coord_xC,lp,cj,sj,Ncj,N)

% Used for symmetric airfoils with zero angle of attack

% Variable Definition
b = zeros(N,1);         % Vector of independent terms
A = zeros(N,N);         % Influence coefficient matrix 
xLoc = zeros(N,N);      % X coordinate of the control point i in the local frame j       
zLoc = zeros(N,N);      % Y coordinate of the control point i in the local frame j
uIndLoc = zeros(N,N);   % Induced horizontal velocity at point i by panel j in the local frame j 
wIndLoc = zeros(N,N);   % Induced vertical velocity at point i by panel j in the local frame j 

uInd = zeros(N,N); % Induced horizontal velocity at point i by panel j in the global frame
wInd = zeros(N,N); % Induced vertical velocity at point i by panel j in the global frame

for i = 1:N
    b(i) = -Ncj(i,:)*Qinf;
    for j = 1:N
        if i == j
            A(i,i) = 0.5;
            uIndLoc(i,i) = 0;
            wIndLoc(i,i) = 0.5;
            uInd(i,i) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
            wInd(i,i) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);
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

            uIndLoc(i,j) = (1/(4*pi))*log((r1^2)/(r2^2));
            wIndLoc(i,j) = (theta2 - theta1)/(2*pi);
            
            % Change of the velocity components to the global frame
            uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
            wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);

            A(i,j) = [uInd(i,j),wInd(i,j)]*Ncj(i,:)';
        end        
    end
end

% Source strenghts obtention 
sigma = A\b;


end