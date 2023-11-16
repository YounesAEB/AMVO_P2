function cp = computeCp(Qinf,V,gamma)
    
    % Control point pressure coefficient computation: Cp_i=1-((abs(Vi)/Qinf)^2
    % cp = 1 - (sqrt(V(:,1).^2+V(:,2).^2)/sqrt(Qinf(1,1)^2+Qinf(2,1)^2)).^2;
    cp = 1 - ((gamma)/sqrt(Qinf(1,1)^2+Qinf(2,1)^2)).^2; 

end