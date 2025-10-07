function [c,f,s]=sav_forest_frontPDEfun(x,t,u,dudx,P)
% Extract parameters

b= P(1);
n=P(2);
m=P(3);
a=P(4);
mu=P(5);

sav=u(1);
fr=u(2);

%PDEs
c=[1;1];
f=[1;0.01].*dudx;
s=[sav*(1-sav)-b*fr*sav-n*sav; mu*fr*(1-fr)-a*sav*fr-m*fr];


