function [cl_int,cl_kutta] = computeCl(cp,lp,Ncj,c,AoA,Qinf,gamma)
    
    % Lift coefficient computation from the pressure coefficient integral:
    cl_int = -1*(sum(cp.*lp.*Ncj)/c) *[-sind(AoA);cosd(AoA)];
    % Lift coefficient computation from Kutta Condition:
    cl_kutta = 2*(sum((gamma.*lp)/(norm(Qinf)*c)));
end