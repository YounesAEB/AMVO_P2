function plotPanelsAndNormVectors(coord_xP,coord_xC,Ncj)
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    figure()
    scatter(coord_xP(:,1),coord_xP(:,2),30,'k','filled');
    hold on
    plot(coord_xP(:,1),coord_xP(:,2),'k');
    quiver(coord_xC(:,1),coord_xC(:,2),Ncj(:,1),Ncj(:,2),'r');
    scatter(coord_xC(:,1),coord_xC(:,2),5,'r');
    hold off

    xlabel('x~(m)')
    ylabel('z~(m)')
    axis equal
    grid on
    grid minor
    axis padded
end