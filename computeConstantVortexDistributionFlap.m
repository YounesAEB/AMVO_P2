function [gamma,uInd,wInd] = computeConstantVortexDistributionFlap(Qinf,coord_xP,coord_xC,lp,cj,sj,Tcj,N,M)

    % Used for asymmetric or symmetric airfoils with non-zero angle of attack
    % Variable Definition
    b = zeros(N+M,1); % Vector of independent terms
    A = zeros(N+M,N+M); % Influence coefficient matrix 
    xLoc = zeros(N+M,N+M); % X coordinate of the control point i in the local frame j       
    zLoc = zeros(N+M,N+M); % Y coordinate of the control point i in the local frame j
    uIndLoc = zeros(N+M,N+M); % Induced horizontal velocity at point i by panel j in the local frame j 
    wIndLoc = zeros(N+M,N+M); % Induced vertical velocity at point i by panel j in the local frame j 
    
    uInd = zeros(N+M,N+M); % Induced horizontal velocity at point i by panel j in the global frame
    wInd = zeros(N+M,N+M); % Induced vertical velocity at point i by panel j in the global frame
    
    for i = 1:N
        b(i) = -Tcj(i,:)*Qinf;
        for j = 1:N
            if i == j
                A(i,i) = -0.5;
                uIndLoc(i,i) = 0.5;
                wIndLoc(i,i) = 0;
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
    
                uIndLoc(i,j) = (theta2 - theta1)/(2*pi);
                wIndLoc(i,j) = (1/(4*pi))*log((r2^2)/(r1^2));
                
                % Change of the velocity components to the global frame
                uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
                wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);
    
                A(i,j) = [uInd(i,j) wInd(i,j)]*Tcj(i,:)';
            end        
        end
        for j = N+1:N+M
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
                wIndLoc(i,j) = (1/(4*pi))*log((r2^2)/(r1^2));
                
                % Change of the velocity components to the global frame
                uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
                wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);
    
                A(i,j) = [uInd(i,j) wInd(i,j)]*Tcj(i,:)';
        end
    end

    for i = N+1:N+M
        b(i) = -Tcj(i,:)*Qinf;
        for j = 1:N
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
                wIndLoc(i,j) = (1/(4*pi))*log((r2^2)/(r1^2));
                
                % Change of the velocity components to the global frame
                uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
                wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);
    
                A(i,j) = [uInd(i,j) wInd(i,j)]*Tcj(i,:)';     
        end
        for j = N+1:N+M
            if i == j
                A(i,i) = -0.5;
                uIndLoc(i,i) = 0.5;
                wIndLoc(i,i) = 0;
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
    
                uIndLoc(i,j) = (theta2 - theta1)/(2*pi);
                wIndLoc(i,j) = (1/(4*pi))*log((r2^2)/(r1^2));
                
                % Change of the velocity components to the global frame
                uInd(i,j) = uIndLoc(i,j)*cj(j) + wIndLoc(i,j)*sj(j);
                wInd(i,j) = -uIndLoc(i,j)*sj(j) + wIndLoc(i,j)*cj(j);
    
                A(i,j) = [uInd(i,j) wInd(i,j)]*Tcj(i,:)';
            end        
        end
    end

    % Kutta condition application (gamma_1 + gamma_N = 0)
    pos1 = round(N/4,0);
    A(pos1,:)   = 0;
    A(pos1,1)   = 1;
    A(pos1,end) = 1;
    b(pos1)     = 0;

    pos2 = N + round(M/4,0);
    A(pos2,:)   = 0;
    A(pos2,1)   = 1;
    A(pos2,end) = 1;
    b(pos2)     = 0;
    
    % Sorce strenghts obtention 
    gamma = A\b;
    gamma(pos1) = (gamma(pos1-1)+gamma(pos1+1))/2;
    gamma(pos2) = (gamma(pos2-1)+gamma(pos2+1))/2;

end