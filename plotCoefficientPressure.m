function plotCoefficientPressure(coord_xC,cp,N,j,AoAaux)
    plot(coord_xC(:,1),cp);
    title("Pressure coefficient plot for " + string(N) + " panels - $C_p$ vs $x/c$")
    xlabel("Chord length $x/c$ ");
    ylabel("Pressure Coefficient $C_{p}$");
    legend("Depends on $\alpha$","Location","southeast");
    %legend("$\alpha=2^\circ$","$\alpha=4^\circ$","$\alpha=6^\circ$","$\alpha=8^\circ$","$\alpha=10^\circ$","Location","southeast");
    grid on;
    grid minor;
    box on;
    axis padded
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);

    if j==numel(AoAaux)
        hold off;
    else
        hold on
    end
end