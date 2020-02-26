function advection2_varv(n)
% example of advection equation in 2D with variable velocity
% n is the number of grid points per dimension (equispaced grid)


if nargin<1, n=64; end
dim=2;
N = [n, n];
o = calc2(N);
X=o.regulargrid;
o.use_twothirds_filtering();


% create initial condition
m = 0*o.gvc(X,1); 

vc{1} =  3*sin(2*o.gvc(X,1)).*cos(2*o.gvc(X,2));
%vc{2} = -cos(2*o.gvc(X,1)).*sin(2*o.gvc(X,2));
vc{2}=vc{1};

% velocity field
v  = o.crv(vc);

% make sure all derivatives return real numbers
o.use_real = true;

flt = o.create_box_filter(1/2);


rhs = @(m) o.fft_filter( - o.inner( o.grad(m), v ) + vc{1}, flt);
%rhs = @(m)  - o.inner( o.grad(m), v ) + vc{1};
%rhs = @(m)-1/2*o.inner(o.grad(m),v)-1/2*o.div( o.scalevec(v,m)) + 1/2*o.div(v).*m + vc{1};

CFL=0.9;
cfl_advect = 1./max(N*norm(o.C(v),inf))
cfltime = CFL*cfl_advect;

monitor_a2([],[],'setup_monitor',[]);
mydata{1}=o;

ops = odeset;
ops = odeset(ops,'RelTol',1e-4,'InitialStep',cfltime);
ops = odeset(ops, 'OutputFcn',@(t,y,flag) monitor_a2(t,y,flag,mydata));
tic;
% o.C() and o.V() are used to convert between vectors and tensor format used in 'calc' class
[T,mt]=ode45(@(t,m) o.C(rhs(o.S(m) )) ,[0,1],o.C(m),ops);
%rk4( @(t,m) o.C(rhs(o.S(m) )),[0,1], o.C(m), ops);
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
