function solution = solve_grda_pde(tend,K)

% solves reaction-diffusion-advection PDE with RHS grda_pde_rhs.m
% and Jacobian Dgrda_pde_rhs.m
% plots solution after each time tend for K iterations

par.Nx=400; % should be even
par.Ny=400; % should be even
%par.Ny=2; % should be even
Nx = par.Nx;
Ny = par.Ny;

%% setup initial condition
%% parameter values from Paul's paper
% par.a=6.2;
% par.b=1;
% par.c = 0.0655;%-0.23;%0.0368;
% par.delta = 100;
% par.m = 1.2;

%% parameter values for circular patch

par.a=6.2;
par.b=1;
par.c =0.0; %0.0655;     %0.115(6.1);%-0.23;%0.0368;
par.delta = 400;
par.m = 1.2;

% par.a=6.2;
% par.b=1;
% par.c =0; %0.0655;     %0.115(6.1);%-0.23;%0.0368;
% par.delta = 100;
% par.m = 1.2;

par.Lx = 400;
Lx = par.Lx;
par.hx = Lx/(Nx-1); hx = par.hx;
x = (0:Nx-1)'*hx;

par.Ly = 400;
Ly = par.Ly;
par.hy = Ly/(Ny-1); hy = par.hy;
y = (0:Ny-1)'*hy;

v1eq=(par.a/par.m+sqrt((par.a/par.m)^2-4*(1+par.a/par.m*par.b)))/(2*(1+par.a/par.m*par.b));
w1eq=par.m*(par.a/par.m-v1eq/(1-par.b*v1eq));
v2eq=0;
w2eq=par.a;

%v0 = repmat([zeros(Nx/2,1);ones(Nx/2,1)*0.55], [Ny,1]);
%v0 = repmat([zeros(Nx/2,1);ones(Nx/2,1)*0.61], [Ny,1]);
%w0 = repmat([par.a*ones(Nx/2,1);ones(Nx/2,1)*2], [Ny,1]);

% v0 = v0+0.02*rand(Nx*Ny,1);
%v0=0.61/2*(1+tanh(x'-Lx/2+2.5*sin(y/4)));
%v0=0.61/2*(1+tanh(x'-Lx/2+2.5*sin(y/4)));

%% intialization with circular patch

%v0=zeros(200);
%w0=zeros(200);

for i=1:400
    for j=1:400
if ((i-200)^2 + (j-200)^2 <= 100^2)
    w0(i,j)=w1eq;  
    v0(i,j)=v1eq;

else
     w0(i,j)=w2eq;  
     v0(i,j)=v2eq;
end
    end
end

v0 = v0+0.02*rand(Nx,Ny);
%%

v0 = reshape(v0',[Nx*Ny,1]);  %% Question: what is the use of this?
w0 = reshape(w0',[Nx*Ny,1]);

sol = [v0;w0];

%% differentiation matrices
e0x = ones(Nx,1);
e0y = ones(Ny,1);

% identity matrices
ex = sparse(1:Nx,[1:Nx],e0x,Nx,Nx); % Nx identity
ey = sparse(1:Ny,[1:Ny],e0y,Ny,Ny); % Ny identity

% d_x
Dx = sparse(1:Nx-1,[2:Nx-1 Nx],ones(Nx-1,1)/2,Nx,Nx);
%Dx(1,Nx) = -1/2; % Periodic boundary conditions
Dx = (Dx - Dx')/hx;
Dx(1,2) = 0; Dx(Nx,Nx-1) = 0; % Neumann boundary conditions

% d_y
Dy = sparse(1:Ny-1,[2:Ny-1 Ny],ones(Ny-1,1)/2,Ny,Ny);
Dy(1,Ny) = -1/2; % Periodic boundary conditions
Dy = (Dy - Dy')/hy;
%Dy(1,2) = 0; Dy(Ny,Ny-1) = 0; % Neumann boundary conditions

% d_xx
D2x = sparse(1:Nx-1,[2:Nx-1 Nx],ones(Nx-1,1),Nx,Nx) - sparse(1:Nx,[1:Nx],e0x,Nx,Nx);
D2x = (D2x + D2x');
%D2x(1,Nx) = 1; D2x(Nx,1) = 1; % Periodic boundary conditions
D2x(1,2)=2; D2x(Nx,Nx-1)=2; % Neumann boundary conditions
D2x = D2x/hx^2;

% d_yy
D2y = sparse(1:Ny-1,[2:Ny-1 Ny],ones(Ny-1,1),Ny,Ny) - sparse(1:Ny,[1:Ny],e0y,Ny,Ny);
D2y = (D2y + D2y');
%D2y(1,2)=2; D2y(Ny,Ny-1)=2; % Neumann boundary conditions
D2y(1,Ny) = 1; D2y(Ny,1) = 1; % Periodic boundary conditions
D2y = D2y/hy^2;

% create differentiation matrices
Dx = sparse(kron(ey,Dx));
Dy = sparse(kron(Dy,ex));
D2x = sparse(kron(ey,D2x));
D2y = sparse(kron(D2y,ex));

par.Dx = Dx;
par.Dy = Dy;
par.D2x = D2x;
par.D2y = D2y;
%par.D2y=0;
%% solve PDE

solution= [];
times = [];
Dgrda_pde_rhs_s = @(t,y)Dgrda_pde_rhs(t,y,par);

%options=odeset('RelTol',1e-8,'AbsTol',1e-8,'Jacobian',Dgrda_pde_rhs_s);
options=odeset('Jacobian',Dgrda_pde_rhs_s);

for j=0:K-1
    j
    time = [0 tend];
    sol = ode15s(@(t,y)grda_pde_rhs(t,y,par), time, sol,options);
    times = [times j*tend];
    sol = sol.y(:,end);
    solution = [solution sol];
end

%% save solution
save(strcat('circularFront_',num2str(tend),'_',num2str(K),'_a_'),'solution');
video=strcat('circularboundary_',num2str(tend),'_',num2str(K));
%% plot solution
 v = VideoWriter(video);
       open(v);
       
for i=1:length(times)
    figure(1)
    %surf(x,y,reshape(solution(1:Nx*Ny,i),[Nx,Ny])')
    imagesc(reshape(solution(1:Nx*Ny,i),[Nx,Ny])');
    set(gca,'YDir','normal');
    colormap(flipud(copper));
    shading flat
    drawnow
    writeVideo(v,getframe(gcf));pause(0.1);
end
close(v);