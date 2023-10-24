function plotPressureCoefficient(coord_xP,coord_xC,Ncj,cp,N)
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    % Pressure Coefficient chord distribution plot
    figure
    plot(coord_xC(:,1),cp);
    grid on
    grid minor
    axis padded
    title('$c_p$ chord distribution plot')
    xlabel('$x_c (m)$')
    ylabel('$c_p$')

    % Pressure Coefficient angle distribution plot
    angle = linspace(0,360,N);
    figure
    plot(angle,cp);
    xlim([0,360]);
    title('$c_p$ angle distribution')
    xlabel('$\alpha$ ($^o$)')
    ylabel('$c_p$')
    grid on
    grid minor
    axis padded

    % Pressure Coefficient airfoil distribution plot
    figure
    for i = 1:N
        Cp_norm = cp.*Ncj;
        if cp(i) < 0
            quiver(coord_xC(i,1),coord_xC(i,2),-Cp_norm(i,1),-Cp_norm(i,2),'r');
            hold on
        else
            quiver(coord_xC(i,1),coord_xC(i,2),Cp_norm(i,1),Cp_norm(i,2),'b');
            hold on 
        end
    end
    plot(coord_xP(:,1),coord_xP(:,2),'k');
    axis equal
    hold off
    title('$c_p$ distribution')
    xlabel('$x$ (m)')
    ylabel('$z$ (m)')
    grid on
    grid minor
    axis padded
end