function u0 = sav_forest_frontICfun(x,P)

n=P(2);
m=P(3);
mu=P(5);


if x <50
    u0 = [1-n;0];
else
    u0 = [0;(mu-m)/mu];

      
end
