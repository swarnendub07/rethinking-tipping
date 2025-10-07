function f = grda_pde_rhs(t,y,par)

% returns right-hand-side of Gilad model

% parameters
Nx = par.Nx;
Ny = par.Ny;
a=par.a;
b=par.b;
c = par.c;
m = par.m;
delta = par.delta;

% solution
v = y(1:Nx*Ny);
w = y(Nx*Ny+1:2*Nx*Ny);

Dx = par.Dx;
D2x = par.D2x;
D2y = par.D2y;

%% RHS Gilad
f1 = (D2x+D2y)*v+c*Dx*v+v.^2.*w.*(1-b*v) - m*v;
f2 = delta*(D2x+D2y)*w+c*Dx*w+a-w-v.^2.*w;

f = [f1; f2];
end


