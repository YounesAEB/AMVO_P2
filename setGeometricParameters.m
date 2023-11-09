function [coord_xP,coord_xC,lp] = setGeometricParameters(c,N,NACA)

    coord_xP = readmatrix("NACA_" + NACA + "_N_" + string(N) + "_coord.txt");
    coord_xP = c*coord_xP(:,2:3);
    coord_xC = zeros(N,2);  % Panel center nodes
    lp = zeros(N,1);        % Panel length definition

    for i = 1:N
        coord_xC(i,1) = (coord_xP(i,1) + coord_xP(i+1,1))/2;
        coord_xC(i,2) = (coord_xP(i,2) + coord_xP(i+1,2))/2;

        lp(i) = sqrt((coord_xP(i+1,1)-coord_xP(i,1))^2 + ...
                     (coord_xP(i+1,2)-coord_xP(i,2))^2);
    end
    
end