%% Start with clean slate
clear all
close all

%% Parameters
% a=1.3;
% b=1.8;
% n=0.4;
% m=0.02;
%a=2.625;
%DuDx=[1;1]
% eps = 0.1;

seed=5000;
rng(seed);

a=5.92;
b=1;
m1=1.2;
delta=200;
saddle=2*m1*(b+sqrt(1+b^2));
%% Uniform steady states

v1eq=(a/m1+sqrt((a/m1)^2-4*(1+a/m1*b)))/(2*(1+a/m1*b));
w1eq=m1*(a/m1-v1eq/(1-b*v1eq));
v2eq=0;
w2eq=a;
unisol=[v1eq,w1eq,v2eq,w2eq];

%% pattern solution




%% Time settings
timepoints =10001;
Endtime = 300000;
tspan = linspace(0,Endtime,timepoints);
t = unique(tspan);

% noisetime=zeros(timepoints);
% for i=1:length(tspan)
%  if rem(i,100)==0
%      noisetime(i)=1;
%  end
% end 
%random=rand(1,length(t)).*(rem(t,50)==0);
gwn=@(t) 0.00005*(rem(t,50)==0);  %0.00005


%% Parameter drift
% mu0 = 0.6;
% mu1 = 0.6;
% 
% mus = @(t) mu0 + (mu1-mu0)*t/Endtime;
% a0 = 5;
% a1 = 2.3;

a0=5.7;
a1=5.94;

as=@(t) a1 -(a1-a0)*t/Endtime;

% dadt=-0.001;
% TimeStep=Endtime/timepoints;

asplot=a1-(a1-a0)*t/Endtime;

%% Spatial settings
L = 1500;
xpoints = 5001;
x = linspace(0,L,xpoints);
m = 0;

%% Initial condition
IC_type = 'turing';

%% Solving
%sol = pdepe(m, @(x,t,u,DuDx) dry_pde(x,t,u,DuDx, as, b, m1,delta), ...
%    @(x) pdeIC(x,L, IC_type,unisol), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t,odeset('Normcontrol', 1, 'Abstol', 10^(-12)));


sol = pdepe(m, @(x,t,u,DuDx) dry_pde(x,t,u,DuDx, as, b, m1,delta,gwn), ...
    @(x) pdeIC(x,L, IC_type,unisol), @(xl,ul,xr,ur,t) pdeBC(xl,ul,xr,ur,t), x, t);


v = sol(:,:,1);
w = sol(:,:,2);
% meanv = sum(v,2)/xpoints; 
% meanw = sum(w(1:length(t),:))/xpoints; 

%% Plotting
%  figure(1)
% surf(x,t,s)
% shading interp
% view(0,90)
% xlabel('$x$', 'Interpreter', 'latex')
% ylabel('$t$', 'Interpreter', 'latex')
% zlabel('$s$', 'Interpreter', 'latex')
% sbar = colorbar();
% ylabel(sbar, '$s$', 'Interpreter', 'latex')
% colormap(flipud(copper))


figure(2)
plot(x(1,:),v(5000,:)) % 5.8610    %5.8250
figure(3)
plot(x(1,:),v(7500,:)) %5.7460     % 5.7625

figure(2)
plot(x(1,:),w(5000,:)) % 5.8610    % 5.8250
figure(3)
plot(x(1,:),w(7500,:)) %5.7460     % 5.7625

% figure(3)
% plot(x(1,:),w(1000,:))
% subplot(2,1,2)
% surf(x,t,f)
% shading interp
% view(0,90)
% xlabel('$x$', 'Interpreter', 'latex')
% ylabel('$t$', 'Interpreter', 'latex')
% zlabel('$f$', 'Interpreter', 'latex')
% fbar = colorbar();
% ylabel(fbar, '$f$', 'Interpreter', 'latex')

% figure(2)
% plot(x,v(timepoints,:))
% figure(3)
% plot(asplot,meanv)

figure(4)
%subplot(2,1,1)
h1=axes
surf(x,asplot,v)
set(h1, 'Ydir', 'reverse')
%surf(x,t,v)
shading interp
view(0,90)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$a$', 'Interpreter', 'latex')
zlabel('$s$', 'Interpreter', 'latex')
axis tight
sbar = colorbar();
ylabel(sbar, '$s$', 'Interpreter', 'latex')
colormap(flipud(copper))



%% Initial condition
function u = pdeIC(x,L, IC_type,unisol)
    v1 = unisol(1);
    v2 = unisol(3);
    w1 = unisol(2);
    w2 = unisol(4); %a

    if strcmp(IC_type, '2-front')
        v = v1 .* ( (x < L/4) + (x > 3*L/4) ) + v2 .* (x > L/4) .* (x < 3*L/4);
        w = w1 .* ( (x < L/4) + (x > 3*L/4) ) + w2 .* ((x > L/4) .* (x < 3*L/4)); 
    elseif strcmp(IC_type, '3-front')
      %  v = v1 .* ( (x < L/5) + (x>4*L/5)) + v2 .* ((x > L/5).* (x < 4*L/5).*((x<2*L/5)+(x>3*L/5)));
        v = v1 .* ( (x < L/5) + (x>4*L/5)+((x>2*L/5).*(x<3*L/5))) + v2 .* ((x > L/5).* (x < 2*L/5))+ v2 .* ((x > 3*L/5).* (x < 4*L/5));
        w = w1 .* ( (x < L/5) + (x>4*L/5)+ ((x>2*L/5).*(x<3*L/5))) + w2 .*((x > L/5).* (x < 2*L/5))+ w2 .* ((x > 3*L/5).*(x < 4*L/5));
    elseif strcmp(IC_type, '1-front')
        v = v1 .* (x < L/2) + v2 .* (x > L/2);
        w = w1 .* (x < L/2) + w2 .* (x > L/2);
%     elseif strcmp(IC_type, 'pat')
%         v=patsolv;
%         w=patsolw;
%                for i=1:1:length(x)
%                 v(i)=patsolv(i);
%                 w(i)=patsolw(i);
%                end
    elseif strcmp(IC_type, 'turing')
%                  init_patv=load('v0_pattern_5_894.mat');
%                  init_patw=load('w0_pattern_5_894.mat');
%                  vp=init_patv.v(end,:);
%                  wp=init_patw.w(end,:);
                 
        for i=1:1:length(x)
                v(i)=v1+0.00001*rand;
                w(i)=w1+0.00001*rand;

%                  v(i)=v1+0.1+0.1*sin(2*pi*22*x/1500);
%                  w(i)=w1+0.1+0.1*sin(2*pi*22*x/1500);
                 
%                  v(i)=vp(i);
%                  w(i)=wp(i);

                
        end

    end
     u = [v;w];
    % save('initial','u');
end

%% PDE function
% PDE has the form
% c u_t = f_x + s
function [c,f,s] = dry_pde(x,t,u,DuDx, as, b, m1,delta,gwn)

    %mu = mus(t) - 0.5 * exp(-x.^2/100); % spatial variation in parameter
    %mu= 0.6-0.5*exp(-(x-50)^2/100);
    %mu= 0.4489+0.05*sin(0.15*x);
    %mu = mus(t) - 0.2 * (abs(x) < 10); % spatial variation (abrupt)
    
    %a=as(t)*exp(-(x-750)^2/500000);

    a=as(t)+0.01*sin(x/250);  %%% hetero simulation to be added in paper
    %a=as(t)*(1+0.000005*x);  

    %a=as(t)*exp(-0.000003*x);

   % a=as(t);
    noise1=gwn(t)*rand;
    noise2=gwn(t)*rand;
 %  [t a]
   % a = a+dadt*TimeStep;
    
    c = [1;1];
    f = [ DuDx(1); delta * DuDx(2)];
    s = [ u(1).^2.*u(2).* (1 - b*u(1)) - m1 * u(1)+noise1*u(1); ...
         a-u(2)-u(1).^2 .* u(2)+noise2*u(2)];
    
end


% f1 = v.^2.*w.*(1-b*v) - m*v;
% f2 = a-w-v.^2.*w;

%% boundary conditions are of the form
% p + q * f = 0
function [pl, ql, pr, qr] = pdeBC(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

