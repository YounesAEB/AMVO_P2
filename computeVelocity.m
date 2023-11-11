function V = computeVelocity(Qinf,gamma,uInd,wInd,N)

% Control point velocity computation: Vi = Qinf + sum_j=1^N(sig_j*V^~_i,j)
    V = zeros(N,2);
    for i=1:N
        sum = zeros(2,1);   
        for j=1:N
            sum = sum + gamma(j)*[uInd(i,j);wInd(i,j)];
        end
        V(i,:) = Qinf + sum;
    end

end