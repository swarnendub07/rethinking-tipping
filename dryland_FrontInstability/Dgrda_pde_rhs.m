function Df = Dgrda_pde_rhs(t,y,par)

% returns Jacobian of right-hand-side of Gilad model

% parameters
Nx = par.Nx;
Ny = par.Ny;

b=par.b;
c = par.c;
delta = par.delta;
m = par.m;

% solution
v = y(1:Nx*Ny);
w = y(Nx*Ny+1:2*Nx*Ny);

Dx = par.Dx;
D2x = par.D2x;
D2y = par.D2y;

Fv = (D2x+D2y)+c*Dx+spdiags(2*v.*w-3*b*w.*v.^2-m,0, Nx*Ny, Nx*Ny);
Fw = spdiags(v.^2.*(1-b*v),0, Nx*Ny, Nx*Ny);
Gv = spdiags(-2*v.*w,0, Nx*Ny, Nx*Ny);
Gw = delta*(D2x+D2y)+c*Dx+spdiags(-1-v.^2,0, Nx*Ny, Nx*Ny);

Df = [Fv Fw; Gv Gw];
end