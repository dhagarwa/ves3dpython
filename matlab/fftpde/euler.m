%%function [omt,o]=euler
clear all; clear globals; clf; clc; clear classes;
addpath ../../lib;
dim=3;
Nx = 128; 
Ny = Nx;
Nz = Nx;
N = [Nx, Ny, Nz];
freq = [1,2,0];
one_over_Re=1e-4;
Tmax = 2;   % time horizon

% SETUP of algebraic and derivative operations for 3D periodic fields using
% equispaced points
o = calc3(N);
o.use_real=true;
%o.use_hou_filtering;
% velocity 
vel = @(om) -o.curl( o.inv_laplacian_vec(om)); 
euler_flux = @(om,u)  o.outer(om,u) - o.outer(u,om);  % flux tensor


%% PROBLEM CASES
%om = calc3_functions(o, 'two_blobs');
%u = calc3_functions(o,'trigonometric_vector'); om=o.curl(u);
%u = calc3_functions(o,'bell_markus'); om = o.curl(u); 
%om = calc3_functions(o,'kerr');  % the classical singularity problem.
om = calc3_functions(o,'two_vortex_tubes'); % from koumousakos-leonard 2011  paper

% exact solutions
exact_solution = 0;
%soltype='trig1'; %trig0, gauss
if exact_solution
    soltype='gauss';
    vals = calc3_functions(o,'euler_exact','solution',soltype);
    om = vals{1};
end

norm(o.C(o.div(om)),inf)
%%
% some of the initial conditions are not divergence free. 
% reproject to ensure divergence free vorticity and velocity.
u = vel(om);
om = o.curl(u);

%%
%clf; vol3d('cdata',o.inner(om,om)), axis tight; axis off; colormap bone; 
%clf; vol3d('cdata',o.div(om)), axis tight; axis off; colormap bone; 
colorbar;
%%

% functions for EDE right hand side
exrhs = @(t) 0;
if exact_solution
   exrhs = @(t) calc3_functions(o,'euler_exact','solution',soltype,'time',t,'rhs_only',1);
end

%Helmholtz projection to set divergence to zero.
divproj =  @(v)  v - o.grad(o.inv_laplacian(o.div(v)));


% EULER equations
%euler_rhs = @(u,om)  - o.hodge_proj(o.div_T(  euler_flux(om, u ) ));
euler_rhs = @(u,om)  - o.div_T(  euler_flux(om, u ) ) -o.scalevec(u,o.div(om));

% NAVIERSTOKES equations
ns_rhs = @(u,om) euler_rhs(u,om) + one_over_Re * o.laplacian_vec(om);
%

CFL=0.9;
cfltime = CFL./max(N*norm(o.C(vel(om)),inf))
monitor([],[],'setup_monitor',[]);
mydata{1}=o; mydata{2}=[];
if exact_solution
  getfirst = @(x) x{1};
  exact_sol =@(time) getfirst(calc3_functions(o,'euler_exact','solution',soltype,'time',time));
  mydata{2}=exact_sol;
end

%%
ops = odeset;
ops = odeset(ops,'RelTol',1e-10,'InitialStep',cfltime);
ops = odeset(ops, 'OutputFcn',@(t,y,flag) monitor(t,y,flag,mydata));
tic;
%[t,omt]=
ode45( @(t,om) o.C( ns_rhs( vel(o.V(om)), o.V(om) ) + exrhs(t)) ,[0,Tmax], o.C(om), ops);
%rk4( @(t,om) o.C( euler_rhs( vel(o.V(om)), o.V(om) ) + exrhs(t)) ,[0,Tmax], o.C(om), ops);
toc
% o.C() and o.V() are used to convert between vectors and tensor format used in calc3

pause;
%% trapezoidal

k = 40;

t0 = t(k);
t1 = t(k+1);
dt = t(k+1)-t(k);
om0 = o.V(omt(k,:));
om1 = o.V(omt(k+1,:));
u0 = vel(om0);

g = om0 + dt/2 * euler_rhs(u0,om0) + dt/2 * (exrhs(t0)+exrhs(t1));
r = om0 - dt/2 * euler_rhs(u0,om0);
r = -r;
mvec = @(xom) xom - dt/2*o.C(( euler_rhs(vel(o.V(xom)),om0)+ euler_rhs(u0,o.V(xom)) ));
p = o.V(gmres(mvec,o.C(g+r),[],1e-8,20));

om1ap = om0+p;

norm(o.C(om1-om1ap),inf)/norm(o.C(om1),inf)


%%
dt = 4*cfltime;
Nt = max(ceil(Tmax/dt),2);

omt_t = zeros([Nt,size(om)]);
omt_t(1,:,:,:,:)=om;
t0 = 0;
u_interp = zeros([2,size(om)]);
ominterp = u_interp;

ops=odeset(ops,'RelTol',1e-6,'InitialStep',dt/5);
newops.useLineSearch=0;
newops.maxLineSearchSteps=2;
newops.linSolver='gmres';
newops.debug=true;
newops.linits = 10;
use_rk_corr=0;
%
monitor(t0,o.C(om),[],mydata);

tic
for j=2:Nt
    t1 = t0+dt;
    om0 = reshape(omt_t(j-1,:,:,:,:),[Nx,Ny,Nz,3]);
    u0 = vel(om0);
    g = om0 + dt/2 * euler_rhs(u0,om0) + dt/2 * (exrhs(t0)+exrhs(t1));
    r = om0 - dt/2 * euler_rhs(u0,om0);
    r = -r;
    mvec = @(xom,om0) xom - dt/2*o.C(( euler_rhs(vel(o.V(xom)),o.V(om0))+ euler_rhs(vel(o.V(om0)),o.V(xom)) ));
%p = o.V(gmres(mvec,o.C(g+r),[],1e-8,10,[],[],o.C(om0)));
%    om1=om0+p;
    
    nwt_res = @(om1) o.C(o.V(om1)-dt/2*euler_rhs(vel(o.V(om1)),o.V(om1)) - g);
    nwt_jac = mvec;
    maxNewtonSteps=1;
    om1 = o.V(smoothNewton(maxNewtonSteps, o.C(om0), nwt_res, nwt_jac,[],newops)); 
    
    for kk=1:1
    if use_rk_corr
        u1 = vel(om1);
        u_interp(1,:,:,:,:)=u0;
        u_interp(2,:,:,:,:)=u1;
        ominterp(1,:,:,:,:)=om0;
        ominterp(2,:,:,:,:)=om1;
        ppu = interp1([t0,t1], u_interp,'spline','pp');
        ppo = interp1([t0,t1], ominterp,'spline','pp');
        funaux =@(om)reshape(om(2,:,:,:,:),[Nx,Ny,Nz,3]);
        ufun=@(t) funaux(ppval(ppu,[t0,t]));
        ofun=@(t) funaux(ppval(ppo,[t0,t]));
        %rk_euler_rhs = @(t,om) o.C(euler_rhs(ufun(t),o.V(om))+exrhs(t));
        residual = @(t) - o.C((om1-om0)/dt - euler_rhs( ufun(t), ofun(t) ) - exrhs(t));
        sdc_euler_rhs = @(t,om) residual(t) +  o.C(euler_rhs(ufun(t),o.V(om)) + 0*euler_rhs(o.V(om),ofun(t)));
        
        disp('BEGIN RK CORRECTION')
        %[tk,omk]=ode45(rk_euler_rhs,[t0,t1], o.C(om0), ops);
        [tk,omk]=ode45(sdc_euler_rhs,[t0,t1], o.C(0*om0));%, ops);
        om1 = om1 +  o.V(omk(end,:));
      end
    end
    
    omt_t(j,:,:,:,:)=om1;
    monitor(t1, o.C(om1),[],mydata);
    t0 = t1;
 end
toc


    
      
        










