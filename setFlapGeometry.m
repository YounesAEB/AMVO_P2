function [coord_xP,coord_xC,lp] = setFlapGeometry(c1,c2,d,df,N,M,NACA)

    R = [cosd(df),-sind(df);sind(df),cosd(df)]; % Rotation Matrix
    [coord_xPN,coord_xCN,lpN] = setGeometricParameters(c1,N,NACA);
    [coord_xPM,coord_xCM,lpM] = setGeometricParameters(c2,M,NACA);
    coord_xCM(:,1) = coord_xCM(:,1) + d;
    coord_xPM(:,1) = coord_xPM(:,1) + d;    
    for i=1:M+1
    coord_xPM(i,:) = coord_xPM(i,:)*R;
    end
    for i=1:M
    coord_xCM(i,:) = coord_xCM(i,:)*R;
    end
    coord_xCM(:,1) = coord_xCM(:,1) + c1;
    coord_xPM(:,1) = coord_xPM(:,1) + c1;

    coord_xP = [coord_xPN;coord_xPM];
    coord_xC = [coord_xCN;coord_xCM];
    lp = [lpN;lpM];
end