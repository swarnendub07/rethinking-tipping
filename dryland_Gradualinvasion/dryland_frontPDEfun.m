function [c,f,s]=dryland_frontPDEfun(x,t,u,dudx,P)
% Extract parameters

b= P(1);
a=P(2);
m1=P(3);

%PDEs
c=[1;1];
f=[1;100].*dudx;
s = [ u(1).^2.*u(2).* (1 - b*u(1)) - m1 * u(1); a-u(2)-u(1).^2 .* u(2)];

