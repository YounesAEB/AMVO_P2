function cm4 = computeCm4(cp,coord_xC,coord_xP,c)
    
    % Momentum coefficient at c/4 computation: 
    cm4 = 0;
    for i=1:size(cp,1)
    cm4 = cm4 + cp(i)*(((coord_xC(i,1)+c/4)*(coord_xP(i+1,1)-coord_xP(i,1)))+(coord_xC(i,2)*(coord_xP(i+1,2)-coord_xP(i,2))))/(c^2);
    end
    %NOTE: it is +c/4 in this case where the cilinder is centered at the
    %origin
end