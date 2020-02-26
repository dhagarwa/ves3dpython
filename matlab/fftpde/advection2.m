function advection2(n)
% example of advection equation in 2D. 
% n is the number of unknowns per dimension.
% in 3D, the only thing that changes is the visualization and setting up N. The rest of the code
% is the same.

if nargin<1, n=64; end
dim=2;

% set lattice dimensions
N = [n, n];

o = calc2(N);

% create points
X=o.regulargrid;

% create initial condition
m = -2*sin(o.gvc(X,1)).*sin(o.gvc(X,2));

% velocity field
v  = o.const_vec([1,1]*2);

% make sure all derivatives return real numbers
o.use_real = true;

% set right hand side for advection diffusion equation
diffusion_coefficient = 1E-9;

rhs = @(m) diffusion_coefficient*o.laplacian(m) - o.inner( o.grad(m), v );

CFL=1;
cfl_advect = 1./max(N*norm(o.C(v),inf))
cfl_diffuse = 1./(max(N)^2)/diffusion_coefficient;
cfltime = CFL*min(cfl_diffuse,cfl_advect);

monitor_a2([],[],'setup_monitor',[]);
mydata{1}=o;

ops = odeset;
ops = odeset(ops,'RelTol',1e-10,'InitialStep',cfltime);
ops = odeset(ops, 'OutputFcn',@(t,y,flag) monitor_a2(t,y,flag,mydata));
tic;
% o.C() and o.V() are used to convert between vectors and tensor format used in 'calc' class
%[T,mt]=ode45(@(t,m) o.C(rhs(o.S(m) )) ,[0,1],o.C(m),ops);
rk4( @(t,m) o.C(rhs(o.S(m) )),[0,1], o.C(m), ops);
toc

%/******************************************************/
function status= monitor_a2(t,m,flag,mydata)
persistent cnt
matfile_stepsize = 100;
pngfile_stepsize = 5;

status = 0;
if strcmp(flag, 'setup_monitor'), cnt = 0; return; end
if strcmp(flag,'done'), return; end;

o=mydata{1};
mi=m;
m = m';
m = o.S(m(1,:));

if mod(cnt, pngfile_stepsize)
    surf(m),axis off, shading interp, view(2), colormap bone;
    pause(0.05);
end
cnt = cnt + 1;
