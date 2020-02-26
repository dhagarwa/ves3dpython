% simple exaple of variable coefficient Laplacian

clear all;
clear globals;
dim = 3;
N = 32;

o = calc3([N,N,N]);
X=o.regulargrid;

contrast = 1e+2;
fr = 1;
K = 1+contrast*abs(sin(fr*o.gvc(X,1))).*abs(sin(fr*o.gvc(X,2))).*abs(sin(fr*o.gvc(X,3)));
cof = mean(K(:));

% variable coefficient laplacian
varlap = @(u) - o.div( o.scalevec( o.grad(u ),K)));

u = 2*sin(o.gvc(X,1)).*cos(o.gvc(X,2)).*sin(o.gvc(X,3));   % exact solution

b= varlap(u);

mvec = @(x)  o.C( varlap( o.S(x)) );
pvec = @(x) -cof*o.C(fftn(o.inv_laplacian(ifftn(o.S(x)))));

TOL = 1E-6;
MAXIT = 10*N;

fprintf('No preconditioner\n');
tic; y1 = gmres(  mvec , b(:), [], TOL, MAXIT); toc
norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

fprintf('Constant coefficient laplace preconditioner\n'); 
tic;y0 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
norm(ifftn(o.S(y0))-u,inf)/norm(u,inf)
