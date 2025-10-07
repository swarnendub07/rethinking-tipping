
clear; close all;

%% Parameters

L=100;
maxt=1000;

P(1)=1.8; %b
P(2)=0.4; %n
P(3)=0.02; %m
P(4)=1.3; %a 
P(5)=0.4; %mu (=0.4489 is the Maxwell point); 

%%

m=0;
t=linspace(0,maxt,800); %tspan
x=linspace(0,L,1500); %xmesh

%%PDEPE
sol=pdepe(m,@sav_forest_frontPDEfun,@sav_forest_frontICfun,@sav_forest_frontBCfun,x,t,[],P);

 u1=sol(:,:,1);
 u2=sol(:,:,2);

 figure(3)
 plot(x,u1(10,:),'color','y', 'linewidth', 2)
 hold on
 plot(x,u2(10,:),'color','g','linewidth', 2)
 ylim([0 1])
 xlim([0 100])

 
xlabel('Distance x','fontsize',20,'fontweight','b','fontname','arial')

a1=area(x,u1(500,:));
a1.FaceColor = 'y';
alpha(a1,.15);

a2=area(x,u2(500,:));
a2.FaceColor = 'g';
alpha(a2,.15);