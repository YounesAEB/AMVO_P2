function cl = computeCl(cp,lp,Ncj,c,AoA)
    
    % Total lift coefficient computation: 
    cl = -1*(sum(cp.*lp.*Ncj)/c)*[-sind(AoA);cosd(AoA)]; 

end