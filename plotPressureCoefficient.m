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
    title("Pressure coefficent distribution over a NACA 0010 at $\alpha=10^\circ$")

    % Pressure Coefficient airfoil distribution plot
    figure
    hold on
    title("Pressure coefficent distribution over a NACA 0010 at $\alpha=10^\circ$")
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
    xlabel("x/c");
ylabel("z/c");
legend("Negative $C_p$","Positive $C_p$","Location","northeast");
grid on;
grid minor;
box on;
axis padded
axis equal
fontsize(13,"points")
hold off;
end