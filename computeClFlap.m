function cl = computeClFlap(cp,lp,Ncj,c,AoA)
    
    % Lift coefficient computation:
    cl = -1*(sum(cp.*lp.*Ncj)/c) *[-sind(AoA);cosd(AoA)];
end