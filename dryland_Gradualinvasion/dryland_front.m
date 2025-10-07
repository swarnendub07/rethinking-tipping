
clear; 
close all;

%% Parameters

L=100;
maxt=2000;

P(1)=1;        %b
P(2)= 6.4;     %a (=6.3513 is the Maxwell point)
P(3)=1.2;      %m

%%
m=0;
t=linspace(0,maxt,800); %tspan
x=linspace(0,L,1000); %xmesh

%%PDEPE
sol=pdepe(m,@dryland_frontPDEfun,@dryland_frontICfun,@dryland_frontBCfun,x,t,[],P);

%sol=xmesh x tspan x variables
u1=sol(:,:,1);
u2=sol(:,:,2);

 
figure(3)
plot(x,u1(10,:),'color','[0.8500 0.3250 0.0980]', 'linewidth', 2)
hold on

xlabel('Distance x','fontsize',20,'fontweight','b','fontname','arial')


figure(3)
a1=area(x,u1(400,:));
a1.FaceColor = '[0.8500 0.5250 0.0980]';
alpha(a1,.15);
