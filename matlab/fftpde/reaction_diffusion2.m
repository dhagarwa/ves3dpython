function diffusion_reaction2(n)
% example of advection equation in 2D. 
% n is the number of unknowns per dimension.

if nargin<1, n=64; end
dim=2;

% set lattice dimensions
N = [n, n];

o = calc2(N);
% create points
X=o.regulargrid;

% create initial condition
% gaussian
sigma = 2*pi*8/n;
c = exp(-1/2 * ((o.gvc(X,1) -pi/4).^2+(o.gvc(X,2)-pi/4).^2)/(sigma.^2) );
c1 = imgaussfilt(10*randn(n,n).^8,2);
c1 = c1/max(c1(:));
%c = c1.*c;

% make sure all derivatives return real numbers
o.use_real = true;

% set right hand side for advection diffusion equation

% material properties
if 0
    contrast=0;
    material=ones(n,n);  
    rows=randi(n,n/4,1);  cols=randi(n,n/4,1); 
    material(rows+1,2:end-1)=contrast; 
    material(2:end-1,cols+1)=contrast; 
else
    xl = [-pi/2.5,0];
    xr = [ pi/2.5,0];
    [tl,rl] = cart2pol(o.gvc(X,1)-xl(1), 2*o.gvc(X,2)-xl(2));
    [tr,rr] = cart2pol(o.gvc(X,1)-xr(1), 2*o.gvc(X,2)-xr(2));
    
    material = ones(n,n);
    material = material.*(rl<pi/3) + material.*(rr<pi/3);
    
    material1 = ones(n,n);
    [tl,rl] = cart2pol(2*o.gvc(X,1)-xl(2), o.gvc(X,2)-xl(1));
    [tr,rr] = cart2pol(2*o.gvc(X,1)-xr(2), o.gvc(X,2)-xr(1));
    material1 = material1.*(rl<pi/4) + material1.*(rr<pi/4);
    
    
    material2 = ones(n,n);
    xr = [0,0];
    [tr,rr] = cart2pol(o.gvc(X,1)-xr(2), o.gvc(X,2)-xr(1));
    material2 = material2.*(rr>0.9*pi);
    material = material+material1 + material2;    
    mysurf(material); axis equal;
end
    
material=double(~material);
save('material','material');

%  smoothing function to avoid alliasing;


o.use_hou_filtering;
%o.use_twothirds_filtering;


%diffusion = @(c) diffcoef * material surfm.* o.laplacian(c);
%diffusion = @(c) diffcoef * material .* o.div( o.scalevec( o.grad(c), material));


smooth = @(u) imgaussfilt(u,1);
%smooth = @(u) u;
material = reshape(max(0,smooth(material)), n,n);
%mysurf(material); 

diffcoef = 1e-3;
diffusion = @(c) diffcoef * o.div( o.scalevec( o.grad(c), material));
rho = 20;
zo = @(c) (min(1,max(0,c))); zo=@(c) max(0,c);
reaction = @(c) rho * material .* c.*(1-c);
%reaction = @(c) smooth(reaction(c));
rhs = @(c) smooth(diffusion(smooth(c))) + (reaction(zo(c))) 
monitor_a2([],[],'setup_monitor',[]);
mydata{1}=o;
ops = odeset;
ops = odeset(ops,'RelTol',1e-3);
ops = odeset(ops, 'OutputFcn',@(t,y,flag) monitor_a2(t,y,flag,mydata));
tic;
[T,mt]=ode45(@(t,c) o.C(rhs(o.S(c) )) ,[0,1],o.C(zo(c.*material)),ops);
fprintf('Max concentration at T is %f\n', max(mt(end,:)))
toc

%/******************************************************/
function status= monitor_a2(t,m,flag,mydata)
persistent cnt
matfile_stepsize = 100;
pngfile_stepsize = 3;

status = 0;
if strcmp(flag, 'setup_monitor'), cnt = 0; return; end
if strcmp(flag,'done'), return; end;

o=mydata{1};
mi=m;
m = m';
m = o.S(m(1,:));

if mod(cnt, pngfile_stepsize)
    surf(m),axis off, shading interp, view(2), colormap bone; axis equal; %set(gca,'zlim',[0 2]); 
    pause(0.1);
end
cnt = cnt + 1;
