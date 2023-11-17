function Mcrit = computeCriticalMach(cp)
    
    Cp_0 = min(cp);
    isentropicExp = 1.4;

    syms Minf2
    eqn_Cp = Cp_0/(sqrt(1-Minf2)+(Cp_0*Minf2)/(2*sqrt(1-Minf2))*(1+Minf2*(isentropicExp-1)/2)) == 2/(isentropicExp*Minf2)*(((2+(isentropicExp-1)*Minf2)/(1+isentropicExp))^(isentropicExp/(isentropicExp-1))-1);
    Mcrit2 = vpasolve(eqn_Cp, Minf2, [0 1]);
    Mcrit = double(sqrt(Mcrit2(1)));

%     % Graphical method
%     auxCP = Cp_0/(sqrt(1-Mcrit^2)+(Cp_0*Mcrit^2)/(2*sqrt(1-Mcrit^2))*(1+Mcrit^2*(isentropicExp-1)/2));
%     vecM = linspace(0,1,100);
%     
%     LaitonesCorrection = Cp_0./(sqrt(1-vecM.^2)+(Cp_0.*vecM.^2)./(2.*sqrt(1-vecM.^2)).*(1+vecM.^2.*(isentropicExp-1)./2));
%     CriticalCp = 2./(isentropicExp.*vecM.^2).*(((2+(isentropicExp-1).*vecM.^2)./(1+isentropicExp)).^(isentropicExp./(isentropicExp-1))-1);
% 
%     figure
%     hold on
%     title("Determination of critical Mach number $M_{cr}$ for $\alpha$= " + string(AoA) + "$^\circ$");
%     %hold on
%     plot(vecM,LaitonesCorrection);
%     plot(vecM,CriticalCp);
%     xlabel("Freestream Mach number $M_{\infty}$");
%     ylabel("Pressure Coefficient $C_{p}$");
%     legend("Laitone's Rule","$C_p^*$=f($M_{cr}$)");
%     grid on;
%     grid minor;
%     box on;
%     axis padded
%     set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
%     xlim([Mcrit-0.2,Mcrit+0.15]);
%     ylim([auxCP-2,auxCP+2]);
%     hold off;


end