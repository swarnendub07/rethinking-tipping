%% Start with clean slate
clear all
close all

%% Parameters
a=6.2;
b=1;
m1=1.2;

%% Uniform steady state

v1eq=(a/m1+sqrt((a/m1)^2-4*(1+a/m1*b)))/(2*(1+a/m1*b));
w1eq=m1*(a/m1-v1eq/(1-b*v1eq));
v2eq=0;
w2eq=a;
unisol=[v1eq,w1eq,v2eq,w2eq];

%% Time settings
timepoints = 3000;
Endtime = 1200;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);

%% Spatial settings
L = 250;
xpoints = 3000;
x = linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = '4-front';

%% Solving


sol = pdepe(m, @(x,t,u,DuDx) dry_pde(x,t,u,DuDx, eps, a, b, m1), ...
    @(x) pdeIC(x,L, IC_type,unisol), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t,odeset('Normcontrol', 'on', 'Abstol', 10^(-12)));
%
v = sol(:,:,1);
w = sol(:,:,2);



%% Plotting
figure()
surf(x,t,v)
shading interp
view(0,90)
xlabel('Space(x)')
ylabel('time (t)')
sbar = colorbar();
ylabel(sbar, 'vegetation')
colormap(flipud(copper))

figure(2)
plot(x,v(10,:),'color','[0.8500 0.3250 0.0980]','linewidth', 2)

hold on

figure(2)
a1=area(x,v(3000,:));
a1.FaceColor = '[0.8500 0.5250 0.0980]';
alpha(a1,.15);
xlabel('Space(x)')
ylabel('vegetation (v)')

figure(3)
plot(x,w(10,:))
hold on

figure(3)
plot(x,w(3000,:))
xlabel('Space(x)')
ylabel('water (w)')


%% Initial condition
function u = pdeIC(x,L, IC_type, unisol)

    v1 = unisol(1); 
    v2 = unisol(3);
    w1 = unisol(2); 
    w2 = unisol(4); 

    if strcmp(IC_type, '2-front')
        
       v = v2 .* ( (x < L/4) + (x > 3*L/4) ) + v1 .* (x > L/4) .* (x < 3*L/4);
       w = w2 .* ( (x < L/4) + (x > 3*L/4) ) + w1 .* ((x > L/4) .* (x < 3*L/4)); 
        
    elseif strcmp(IC_type, '4-front')
     
        v = v2 .* ( (x < L/5) + (x>4*L/5)+((x>2*L/5).*(x<3*L/5))) + v1 .* ((x > L/5).* (x < 2*L/5))+ v1 .* ((x > 3*L/5).* (x < 4*L/5));
        w = w2 .* ( (x < L/5) + (x>4*L/5)+ ((x>2*L/5).*(x<3*L/5))) + w1 .*((x > L/5).* (x < 2*L/5))+ w1 .* ((x > 3*L/5).*(x < 4*L/5));

    elseif strcmp(IC_type, '1-front')
        v = v1 .* (x < L/2) + v2 .* (x > L/2);
        w = w1 .* (x < L/2) + w2 .* (x > L/2);
    end

    u = [v;w];
end

%% PDE function

function [c,f,s] = dry_pde(x,t,u,DuDx, eps, a, b, m1)


    c = [1;1];
    f = [ DuDx(1); 100 * DuDx(2)];                           % delta_w=100;
    s = [ u(1).^2.*u(2).* (1 - b*u(1)) - m1 * u(1); ...
         a-u(2) - u(1).^2 .* u(2)];
end


%% boundary conditions 

function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

