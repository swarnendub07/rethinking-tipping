function u0 =dryland_frontICfun(x,P)

a=P(2);
b=P(1);
m1=P(3);

v1eq=(a/m1+sqrt((a/m1)^2-4*(1+a/m1*b)))/(2*(1+a/m1*b));
w1eq=m1*(a/m1-v1eq/(1-b*v1eq));

if x <50
    u0 = [0;a];
else
    u0 = [v1eq;w1eq];
      
end
