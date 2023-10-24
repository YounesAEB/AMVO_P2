function plotVelocityDistribution(Qinf,V,N)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    % Velocity distribution plot
    angle = linspace(0,360,N);
    figure
    plot(angle,(V(:,1)/norm(Qinf)),'b');
    hold on
    plot(angle,(V(:,2)/norm(Qinf)),'r');
    hold off
    xlim([0,360]);
    xlabel('$\theta (^o)$')
    ylabel('$V (m/s¿?)$')
    legend('$V_x$','$V_z$',Location='southwest')

    grid on
    grid minor
    axis padded
end