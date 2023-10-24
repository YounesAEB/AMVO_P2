function [cj,sj,Ncj,Tcj] = computePanelAngleAndNormalAndTangentVectors(coord_xP,lp,N)

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

end