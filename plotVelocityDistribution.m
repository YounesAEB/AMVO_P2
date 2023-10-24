function plotVelocityDistribution(Qinf,V,N)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    % Velocity distribution plot
    angle = linspace(0,360,N);
    figure
    plot(angle,(V(:,1)/norm(Qinf)),'r');
    hold on
    plot(angle,(V(:,2)/norm(Qinf)),'b');
    hold off
    xlim([0,360]);

    grid on
    grid minor
    axis padded
end