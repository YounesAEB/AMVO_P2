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

% V = zeros(N,2);
% for i=1:N
% sum_x = 0;
% sum_z = 0;
%     for j=1:N
%         sum_x = sum_x + gamma(j)*uInd(i,j);
%         sum_z = sum_z + gamma(j)*wInd(i,j);
%     end
% V(i,1) = Qinf(1,1) + sum_x;
% V(i,2) = Qinf(2,1) + sum_z;
% end

end