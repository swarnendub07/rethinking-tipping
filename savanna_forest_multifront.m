%% Start with clean slate
clear all
close all

%% Parameters
a=1.3;
b=1.8;
n=0.4;
m1=0.02;
mu=0.46;
eps = 0.1;

%% Uniform steady state

s1eq=0;
f1eq=(mu-m1)/mu;
s2eq=1-n;
f2eq=0;
unisol=[s1eq,f1eq,s2eq,f2eq];

%% Time settings
timepoints = 1500;
Endtime = 14000;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);

%% Parameter drift
mu0 = 0.46;
mu1 = 0.46;
mus = @(t) mu0 + (mu1-mu0)*t/Endtime;


%% Spatial settings
L = 250;
xpoints = 1000;
x = linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = '2-front';

%% Solving
sol = pdepe(m, @(x,t,u,DuDx) humid_pde(x,t,u,DuDx, eps, a, b, n, m1, mus), ...
    @(x) pdeIC(x,L, IC_type,unisol), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t);

s = sol(:,:,1);
f = sol(:,:,2);

%% Plotting
figure(1)
surf(x,t,s)
shading interp
view(0,90)
xlabel('Space(x)')
ylabel('time(t)')
sbar = colorbar();
ylabel(sbar, 'savanna')
colormap(summer)

figure(2)
surf(x,t,f)
shading interp
view(0,90)
xlabel('Space(x)')
ylabel('time(t)')
fbar = colorbar();
ylabel(fbar, 'forest')
colormap(flipud(summer))


figure(3)
a1=area(x,s(1500,:));
a1.FaceColor = 'y';
alpha(a1,.15);

hold on

figure(3)
a2=area(x,f(1500,:));
a2.FaceColor = 'g';
alpha(a2,.15);
xlabel('Space(x)')
ylabel('savanna/forest')

%% Initial condition
function u = pdeIC(x,L, IC_type,unisol)
    % s1 = 0;
    % s2 = 0.6;
    % f1 = 0.9554;
    % f2 = 0;

    s1 = unisol(1); 
    s2 = unisol(3);
    f1 = unisol(2); 
    f2 = unisol(4); 

    if strcmp(IC_type, '2-front')
        s = s1 .* ( (x < L/4) + (x > 3*L/4) ) + s2 .* (x > L/4) .* (x < 3*L/4);
        f = f1 .* ( (x < L/4) + (x > 3*L/4) ) + f2 .* ((x > L/4) .* (x < 3*L/4)); 
    elseif strcmp(IC_type, '4-front')
        s = s1 .* ( (x < L/5) + (x>4*L/5)) + s2 .* ((x > L/5).* (x < 4*L/5).*((x<2*L/5)+(x>3*L/5)));
        f = f1 .* ( (x < L/5) + (x>4*L/5)+ ((x>2*L/5).*(x<3*L/5))) + f2 .* ((x > L/5).* (x < 4*L/5));
    elseif strcmp(IC_type, '1-front')
        s = s1 .* (x < L/2) + s2 .* (x > L/2);
        f = f1 .* (x < L/2) + f2 .* (x > L/2);
    end

    u = [s;f];
end

%% PDE function

function [c,f,s] = humid_pde(x,t,u,DuDx, eps, a, b, n, m1,mus)
    
    mu= mus(t)+0.025*sin(0.16*x);  % with spatial heterogeneity
       
    % mu=mus(t);                     % without spatial heterogeneity
    
    c = [1;1];
    f = [ DuDx(1); eps.^2 * DuDx(2)];                           % eps^2=delta=0.01;
    s = [ u(1) .* (1 - u(1)) - b .* u(1) .* u(2) - n * u(1); ...
        mu .* u(2) .* (1 - u(2)) - a .* u(1) .* u(2) - m1 .* u(2)];
end

%% boundary conditions

function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

