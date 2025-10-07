function[pl,ql,pr,qr] = sav_forest_frontBCfun(xl,ul,xr,ur,t,P)
% No flux BCs on both sides
pl = [0;0]; ql = [1;1];
pr = [0;0]; qr = [1;1];
