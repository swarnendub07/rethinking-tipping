
%% Start with clean slate
clear all
close all

load("v_init.mat");
load("w_init.mat");


patsol=[v; w];

%% Parameter for the dryland model

b=1;
m1=0.9; 
a=4.6;

%% Uniform steady state

v1eq=(a/m1+sqrt((a/m1)^2-4*(1+a/m1*b)))/(2*(1+a/m1*b));
w1eq=m1*(a/m1-v1eq/(1-b*v1eq));
v2eq=0;
w2eq=a;
unisol=[v1eq,w1eq,v2eq,w2eq];

%% Time settings
timepoints = 10000;
Endtime = 10000;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);

%% Parameter drift

a0=4.6;
a1=4.15;

as=@(t) a0 -(a0-a1)*t/Endtime;
asplot=(a0-(a0-a1)*t/Endtime);


%% Spatial settings
L = 200;
xpoints = 3000;
x = linspace(0,L,xpoints);
xvec=x;
%x1= linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = '4-front';

%% Solving
sol = pdepe(m, @(x,t,u,DuDx) dry_pde(x,t,u,DuDx, as, b, m1), ...
    @(x) pdeIC(x,L, IC_type,unisol, patsol,xvec), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t,odeset('Normcontrol', 'on', 'AbsTol', 10^(-12)));
%
v = sol(:,:,1);
w = sol(:,:,2);

%% Plotting
figure()
surf(x,asplot,v)
%surf(x,t,v)
shading interp
view(0,90)
xlabel('Space(x)')
ylabel('time(t)')
sbar = colorbar();
ylabel(sbar, 'vegetation')
colormap(flipud(copper))
set(gca, 'YDir','reverse')

figure(2)
plot(x,v(1,:),'color','[0.8500 0.3250 0.0980]','linewidth', 2)
hold on

figure(2)
a2=area(x,v(8000,:));
a2.FaceColor = '[0.8500 0.5250 0.0980]';
alpha(a2,.15);


figure(3)
plot(x,w(10,:))
hold on
figure(3)
plot(x,w(8000,:))


%% Initial condition
function u = pdeIC(x,L, IC_type,unisol,patsol,xvec)

    v1 = unisol(1); 
    v2 = unisol(3);
    w1 = unisol(2); 
    w2 = unisol(4); 
     
    if strcmp(IC_type, '2-front')
        
        v = v2 .* ( (x < L/4) + (x > 3*L/4) ) + v1 .* (x > L/4) .* (x < 3*L/4);
        w = w2 .* ( (x < L/4) + (x > 3*L/4) ) + w1 .* ((x > L/4) .* (x < 3*L/4)); 
       
    elseif strcmp(IC_type, '4-front')

         [R,C]=find(xvec==x,1);
         v=patsol(1,C);
         w=patsol(2,C);
         
         % v = v2 .* ( (x < L/5) + (x>4*L/5)+((x>2*L/5).*(x<3*L/5))) + v1 .* ((x > L/5).* (x < 2*L/5))+ v1 .* ((x > 3*L/5).* (x < 4*L/5));
         % w = w2 .* ( (x < L/5) + (x>4*L/5)+ ((x>2*L/5).*(x<3*L/5))) + w1 .*((x > L/5).* (x < 2*L/5))+ w1 .* ((x > 3*L/5).*(x < 4*L/5));
              
         elseif strcmp(IC_type, '1-front')
         v = v1 .* (x < L/2) + v2 .* (x > L/2);
         w = w1 .* (x < L/2) + w2 .* (x > L/2);  
         
    end

    u = [v;w];
end

%% PDE function


function [c,f,s] = dry_pde(x,t,u,DuDx, as, b, m1)

     a=as(t);

    c = [1;1];
    f = [ DuDx(1); 2500 * DuDx(2)];                    % delta_w=2500;
    s = [ u(1).^2.*u(2).* (1 - b*u(1)) - m1 * u(1); ...
         a-u(2) - u(1).^2 .* u(2)];
end   


%% boundary conditions are of the form

function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

