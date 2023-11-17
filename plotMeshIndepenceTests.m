function plotMeshIndepenceTests(AoAaux,Naux,CL,CM4,MCR)
    
    figure
    hold on
    title("Panel density independence test - $C_l$ vs $\alpha$")
    for i=1:numel(Naux)
        plot(AoAaux,CL(i,:));
    end
    xlabel("Angle of attack $\alpha$ ($^\circ$)");
    ylabel("Lift Coefficient $C_{l}$");
    legend("$Depends on Nelem$","Location","northwest");
    %legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","northwest");
    grid on;
    grid minor;
    box on;
    axis padded
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
    hold off;
    
    figure
    hold on
    title("Panel density independence test - ${C_m}_{1/4}$ vs $\alpha$")
    for i=1:numel(Naux)
        plot(AoAaux,CM4(i,:));
    end
    xlabel("Angle of attack $\alpha$ ($^\circ$)");
    ylabel("Momentum Coefficient about the quarter ${C_{m}}_{1/4}$");
    legend("$Depends on Nelem$","Location","southwest");
    %legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","southwest");
    grid on;
    grid minor;
    box on;
    axis padded
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
    hold off;
    
    figure
    hold on
    title("Panel density independence test - $M_{crit}$ vs $\alpha$")
    for i=1:numel(Naux)
        plot(AoAaux,MCR(i,:));
    end
    xlabel("Angle of attack $\alpha$ ($^\circ$)");
    ylabel("Critical Freestream Mach number $M_{crit}$");
    legend("$Depends on Nelem$","Location","northwest");
    %legend("$N_{elem}=16$","$N_{elem}=32$","$N_{elem}=64$","$N_{elem}=128$","$N_{elem}=256$","$N_{elem}=512$","Location","southwest");
    grid on;
    grid minor;
    box on;
    axis padded
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
    hold off;
end