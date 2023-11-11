function cm4 = computeCm4Flap(cp,coord_xC,coord_xP,c,N,M)
    
    % Momentum coefficient at c/4 computation: 
    cm4 = 0;
    for i=1:N
    cm4 = cm4 + cp(i)*(((coord_xC(i,1)-c/4)*(coord_xP(i+1,1)-coord_xP(i,1)))+(coord_xC(i,2)*(coord_xP(i+1,2)-coord_xP(i,2))))/(c^2);
    end
    for i=N+1:N+M
    cm4 = cm4 + cp(i)*(((coord_xC(i,1)-c/4)*(coord_xP(i+2,1)-coord_xP(i+1,1)))+(coord_xC(i,2)*(coord_xP(i+2,2)-coord_xP(i+1,2))))/(c^2);
    end
end