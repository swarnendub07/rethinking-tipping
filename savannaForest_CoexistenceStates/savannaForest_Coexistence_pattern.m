%% Start with clean slate

clear all
close all

%% Parameters
a=1.3;
b=1.8;
n=0.4;
m1=0.02;                     % m
eps = 0.01;                  % \delta
mu_maxwell=0.4489;           % Maxwell point (\mu_M)

%% Time settings
timepoints = 1500;
Endtime = 14000;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);

%% Forest tree growth rate
mus=0.46;

%% Spatial settings
L = 250;
xpoints = 1000;
x = linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = '2-front';

%% Characteristic of the spatial heterogeneity

amp= 0.025;    
wavenum=0.05;  % wavenum=0.16 has been used in the Supplemental PDF

%% Determines location of the initial savanna domain

sizeLeft=0;    %Fig 6B :sizeLeft=45, sizeRight=5;
sizeRight=0;   %Fig 6C :sizeLeft=45, sizRight=-15;
               %Fig 6D :sizeLeft=60, sizRight=0;

%% solving

sol = pdepe(m, @(x,t,u,DuDx) humid_pde(x,t,u,DuDx, eps, a, b, n, m1, amp, wavenum, mus), ...
    @(x) pdeIC(x,L, IC_type,sizeLeft,sizeRight), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t);

s = sol(:,:,1);
f = sol(:,:,2);

%% Plotting the savanna with time

figure(1)
surf(x,t,s)
shading interp
view(0,90)
xlabel('Space(x)','FontSize',18, FontWeight='bold')
ylabel('time(t)','FontSize',18, FontWeight='bold')
zlabel('$s$', 'Interpreter', 'latex')
sbar = colorbar();
ylabel(sbar, 'savanna', 'FontSize',18, FontWeight='bold')
colormap(summer)

%% Plotting the forest with time
figure(2)
surf(x,t,f)
shading interp
view(0,90)
xlabel('Space(x)','FontSize',18, FontWeight='bold')
ylabel('time(t)','FontSize',18, FontWeight='bold')
zlabel('$f$', 'Interpreter', 'latex')
xlim([0 250])
ylim([0 14000])
yticks([0 14000])
xticks([0 250])
fbar = colorbar();
ylabel(fbar, 'forest','FontSize',18, FontWeight='bold')
colormap(flipud(summer))
ax = gca;          
ax.FontSize = 18; 
ax.FontWeight = 'bold'; 

hold on 

%% calculating the locations favourable for forest or savanna invasion
%% \mu =0.46;

mu_vals = 0.46 + amp * sin(wavenum * x);
diff=mu_vals-mu_maxwell;
tolerance = 0.5e-3;
indices = find(abs(mu_vals - mu_maxwell) < tolerance);
x_solution = x(indices);

%% plotting the locations favourable for forest or savanna invasion

yyaxis left
diff_pos = diff;
diff_pos(diff_pos <= 0) = NaN;

yyaxis right
plot(x, diff_pos, LineWidth=2, Color=[0, 0.392, 0]);
ylim([0 max(diff_pos) * 10]);

ax1 = gca; 
ax1.FontSize = 16; 
ax1.FontName = 'Arial';
ax1.FontWeight = 'bold'; 
ax1.YTick = [];

ax2 = axes('Position', ax1.Position, 'Color', 'none', ...
           'YAxisLocation', 'right', 'XTick', [], 'XColor', 'k', 'YColor', 'k');
hold(ax2, 'on');
ax2.FontSize = 16; 
ax2.FontWeight = 'bold'; 
ax2.FontName = 'Arial';

diff_neg = diff;
diff_neg(diff_neg >= 0) = NaN;

plot(ax2, x, diff_neg, 'LineWidth', 2, Color=[0.8, 0.6, 0]); 
ax2.YLim = [min(diff_neg)*10 0];  % min(diff_neg) is negative
ax2.YTick = [];
ax2.XLim = ax1.XLim;   % sync x range
linkaxes([ax1 ax2], 'x'); % keep them linked  
linkprop([ax1 ax2], 'Position');  % synchronize x-axes

% Set figure size
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 12, 8.5])
%exportgraphics(gcf, 'no_hetero_sL0_sR0.tif', 'Resolution', 300)


%% plotting the savanna and forest density at the last time point

figure(3)
a1=area(x,s(1500,:));
a1.FaceColor = 'y';
alpha(a1,.15);
hold on

figure(3)
a2=area(x,f(1500,:));
a2.FaceColor = 'g';
alpha(a2,.15);
hold on


%% Plotting the spatial heterogeneity
figure(4)
plot(x, mu_vals,  'LineWidth', 2,'Color','b');
line([0 250], [0.46 0.46], 'Color','b','LineWidth', 2, LineStyle='--')
hold on
line([0 250], [mu_maxwell mu_maxwell], 'Color','k','LineWidth', 2)
ylim ([0.43 0.49])
yticks([0.43 0.49])
xlim([0 250])
xticks([0 250])
ax = gca;          % get current axes
ax.FontSize = 20; 
ax.FontName = 'Arial'; 
ax.FontWeight = 'bold'; 
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 12, 6])
%exportgraphics(gcf, 'spatial_hetero_wav_0_05.tif', 'Resolution', 300)
%savefig(gcf, 'spatial_hetero_wav_0_05.fig', 'Resolution', 300)

%% Initial condition
function u = pdeIC(x,L, IC_type,sizeLeft,sizeRight)
    s1 = 0;
    s2 = 0.6;
    f1 = 0.9554;
    f2 = 0;

    if strcmp(IC_type, '2-front')
        s = s1 .* ( (x < L/4+sizeLeft) + (x > 3*L/4-sizeRight) ) + s2 .* (x > L/4+sizeLeft) .* (x < 3*L/4-sizeRight);
        f = f1 .* ( (x < L/4+sizeLeft) + (x > 3*L/4-sizeRight) ) + f2 .* ((x > L/4+sizeLeft) .* (x < 3*L/4-sizeRight)); 
    elseif strcmp(IC_type, '3-front')
        s = s1 .* ( (x < L/5) + (x>4*L/5)) + s2 .* ((x > L/5).* (x < 4*L/5).*((x<2*L/5)+(x>3*L/5)));
        f = f1 .* ( (x < L/5) + (x>4*L/5)+ ((x>2*L/5).*(x<3*L/5))) + f2 .* ((x > L/5).* (x < 4*L/5));
    elseif strcmp(IC_type, '1-front')
        s = s1 .* (x < L/2) + s2 .* (x > L/2);
        f = f1 .* (x < L/2) + f2 .* (x > L/2);        
    end

    u = [s;f];
end

%% PDE function
% PDE has the form
% c u_t = f_x + s
function [c,f,s] = humid_pde(x,t,u,DuDx, eps, a, b, n, m1,amp,wavenum,mus)

    
     mu= mus+amp*sin(wavenum*x); % spatial environmental heterogeneity
 %   mu=mus;                     % absence of spatial environmental heterogeneity
    
   
    c = [1;1];
    f = [ DuDx(1); eps * DuDx(2)];
    s = [ u(1) .* (1 - u(1)) - b .* u(1) .* u(2) - n * u(1); ...
        mu .* u(2) .* (1 - u(2)) - a .* u(1) .* u(2) - m1 .* u(2)];
end

%% boundary conditions are of the form
% p + q * f = 0
function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

