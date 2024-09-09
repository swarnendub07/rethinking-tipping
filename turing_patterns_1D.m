%% Start with clean slate
clear all
close all

%% Parameter for the dryland model

a=5.95;
b=1;
m1=1.2;
delta=200;

%% tipping point in the non-spatial dryland model

%saddle=2*m1*(b+sqrt(1+b^2));


%% Uniform steady states

v1eq=(a/m1+sqrt((a/m1)^2-4*(1+a/m1*b)))/(2*(1+a/m1*b));
w1eq=m1*(a/m1-v1eq/(1-b*v1eq));
v2eq=0;
w2eq=a;
unisol=[v1eq,w1eq,v2eq,w2eq];

%% Time settings
timepoints =10001;
Endtime = 200000;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);


gwn=@(t) 0.00005*(rem(t,50)==0);    % creating a vector with small intensity of noise
                                    % after every 50 time points


%% Parameter drift for dryland

a0=5.7;
a1=5.95;
as=@(t) a1 -(a1-a0)*t/Endtime;

asplot=a1-(a1-a0)*t/Endtime;        % for plotting 'a' as one of the axes

%% Spatial settings
L = 1500;
xpoints = 5001;
x = linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = 'turing';           

%% Solving the pde

sol = pdepe(m, @(x,t,u,DuDx) dry_pde(x,t,u,DuDx, as, b, m1,delta,gwn), ...
    @(x) pdeIC(x,L, IC_type,unisol), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t);


v = sol(:,:,1);                     % vegetation
w = sol(:,:,2);                     % water



%% Plotting

figure(1)
plot(x(1,:),v(5000,:))
xlabel('Space(x)')
ylabel('vegetation')
figure(2)
plot(x(1,:),v(7500,:))
xlabel('Space(x)')
ylabel('vegetation')


figure(3)
plot(x(1,:),w(5000,:))
xlabel('Space(x)')
ylabel('water')
figure(4)
plot(x(1,:),w(7500,:)) 
xlabel('Space(x)')
ylabel('water')


figure(5)
h1=axes;
surf(x,asplot,v)
set(h1, 'Ydir', 'reverse')
%surf(x,t,v)
shading interp
view(0,90)
xlabel('Space(x)')
ylabel('p')
axis tight
sbar = colorbar();
ylabel(sbar, 'vegetation')
colormap(flipud(copper))



%% Initial condition
function u = pdeIC(x,L, IC_type,unisol)
    v1 = unisol(1);
    v2 = unisol(3);
    w1 = unisol(2);
    w2 = unisol(4);

    if strcmp(IC_type, '2-front')
        v = v1 .* ( (x < L/4) + (x > 3*L/4) ) + v2 .* (x > L/4) .* (x < 3*L/4);
        w = w1 .* ( (x < L/4) + (x > 3*L/4) ) + w2 .* ((x > L/4) .* (x < 3*L/4)); 
    elseif strcmp(IC_type, '1-front')
        v = v1 .* (x < L/2) + v2 .* (x > L/2);
        w = w1 .* (x < L/2) + w2 .* (x > L/2);
    elseif strcmp(IC_type, 'turing')
                 
        for i=1:1:length(x)
                v(i)=v1+0.00001*rand;              % adding small perturbation to the uniform steady state
                w(i)=w1+0.00001*rand;
                
        end

    end
     u = [v;w];
end

%% PDE function

function [c,f,s] = dry_pde(x,t,u,DuDx, as, b, m1,delta,gwn)


    a=as(t);
    noise1=gwn(t)*rand;
    noise2=gwn(t)*rand;

    c = [1;1];

    f = [ DuDx(1); delta * DuDx(2)];
    s = [ u(1).^2.*u(2).* (1 - b*u(1)) - m1 * u(1)+noise1*u(1); ...
         a-u(2)-u(1).^2 .* u(2)+noise2*u(2)];
    
end

%% boundary conditions

function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

