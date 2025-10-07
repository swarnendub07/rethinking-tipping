clear all
clc

%% When \mu=0.4 (Supplemental PDF)

%eps=[0.0025 0.005 0.01 0.05];           
%Maxwell=[0.4662 0.4591 0.4489 0.40595]; 

%% When \mu=0.46 (Main text)

eps=[0.005 0.01 0.05 1 5 25];     %\delta
Maxwell=[0.4591 0.4489 0.40595 0.25824 0.19355 0.16053]; % Maxwell point (\mu_M)

%% plotting the Maxwell point (\mu_M) with respect to \delta

yyaxis left
semilogx(eps, Maxwell, '-o', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6, 'Color','k', LineWidth=2);
ylabel('Maxwell point (\mu_M)', 'FontSize',18, FontWeight='bold', Color='k');
%ylim([0.35,0.5])
ylim([0.1,0.5])
yticks([0.1, 0.2, 0.3, 0.4, 0.46, 0.5])
yticklabels({'0.1','0.2','0.3','0.4','\mu','0.5'})
ax1 = gca;
ax1.YColor = 'k'; 

hold on

%line([0.001 100],[0.4,0.4], 'Color','b', LineWidth=2)
line([0.001 100],[0.46,0.46], 'Color','b', LineWidth=2)

%% plotting the Threshold amplitude required for coexistence state 
%% (\mu_Ac=|Maxwell - 0.46|) with respect to \delta
 
yyaxis right
%semilogx(eps, abs(Maxwell - 0.4), '-o', ...
    %'MarkerFaceColor', 'r', 'MarkerSize', 6, 'Color','r', LineWidth=2);
semilogx(eps, abs(Maxwell - 0.46), '-o', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6, 'Color','r', LineWidth=2);
ylabel('Threshold ampl. (\mu_{Ac})', FontSize=18, FontWeight='bold',Color='r');
%ylim([0 0.1])
ylim([0 0.35])

xlabel(sprintf('Ratio of diffusion coefficients\n(\\delta)'), ...
    'FontSize',18, 'FontWeight','bold');
%xlim([10^(-3) 10^(-1)])
xlim([10^(-3) 10^(2)])

%xticks([10^(-3) 10^(-2) 10^(-1)])
xticks([10^(-3) 10^(-2) 10^(-1) 1 10^(1) 10^(2)])


ax=gca;
ax.YColor = 'r'; 
ax.FontSize = 18; 
ax.FontWeight = 'bold'; 

set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 12, 10])


%exportgraphics(gcf, 'Threshold_amplitude_0_46.tif', 'Resolution', 300)
