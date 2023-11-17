function plotExperimentalAndPanelMethodCl(AoAaux,CL_kutta)
    
    NACA0010_Re100k_data = readmatrix("xf-naca0010-il-100000-n5.txt");
    NACA0010_Re500k_data = readmatrix("xf-naca0010-il-500000-n5.txt");
    NACA0010_Re1M_data = readmatrix("xf-naca0010-il-1000000-n5.txt");
    
    AoAExp_Re100k = NACA0010_Re100k_data(41:89,1);
    ClExp_Re100k = NACA0010_Re100k_data(41:89,2);
    
    AoAExp_Re500k = NACA0010_Re500k_data(:,1);
    ClExp_Re500k = NACA0010_Re500k_data(:,2);
    
    AoAExp_Re1M = NACA0010_Re1M_data(:,1);
    ClExp_Re1M = NACA0010_Re1M_data(:,2);

    figure
    hold on
    title("DVM vs experimental data")
    scatter(AoAaux,CL_kutta)
    %plot(AoAExp_Re100k,ClExp_Re100k);
    plot(AoAExp_Re500k,ClExp_Re500k);
    plot(AoAExp_Re1M,ClExp_Re1M);
    xlabel("Angle of attack $\alpha$ ($^\circ$)");
    ylabel("Lift Coefficient $C_{l}$");
    legend("DVM","Experimental $Re=5\cdot 10^5$", "Experimental $Re=1\cdot 10^6$","Location","southeast");
    grid on;
    grid minor;
    box on;
    axis padded
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
    hold off;


%     % Exportar a PDF 
%     set(gcf, 'Units', 'Centimeters');
%     pos = get(gcf, 'Position');
%     set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
%         'PaperSize',[pos(3), pos(4)]);
%     print(gcf, 'ExpPMCl', '-dpdf', '-r0'); % incrementar '-r0' resolución 
end

