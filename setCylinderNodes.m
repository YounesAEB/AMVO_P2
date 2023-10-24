function [coord_xP] = setCylinderNodes(R,N)
    
    coord_xP = zeros(N+1,2); % Panel limits coordinates
    
    for i = 1:N+1
        coord_xP(i,1) = R*cosd(-(360/N)*(i-1));
        coord_xP(i,2) = R*sind(-(360/N)*(i-1));
    end
    
end