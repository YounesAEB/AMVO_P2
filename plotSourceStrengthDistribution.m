function plotSourceStrengthDistribution(coord_xC,coord_xP,Ncj,gamma,N)

    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    intensity = zeros(N,2);
    f         = 0.05; %Scale Factor
    figure
    for i = 1:N
        intensity(i,:) = f*gamma(i)*Ncj(i,:);
        if gamma(i) < 0
            quiver(coord_xC(i,1),coord_xC(i,2),-intensity(i,1),-intensity(i,2),'r');
            hold on
        else
            quiver(coord_xC(i,1),coord_xC(i,2),intensity(i,1),intensity(i,2),'b');
            hold on 
        end
    end
    plot(coord_xP(:,1),coord_xP(:,2),'k');
    hold off
    
    axis equal
    grid on
    grid minor
    axis padded
end