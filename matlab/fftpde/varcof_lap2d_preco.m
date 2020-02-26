% study preconditioners (laplacian, sparsified, simple HSS) for variable coefficient laplacian.

clear all;
clear globals;
dim = 2;
N = 32;

o = calc2([N,N]);
id = eye(N*N);
X=o.regulargrid;

u = 2*sin(o.gvc(X,1)).*cos(o.gvc(X,2));

% It is better to have the Laplacian in the frequency domain for sparsity
lap = @(u)   - fftn(o.div( o.grad(ifftn(u))));

contrast = 1e+4;
fr = 1;
K = 1+contrast*abs(sin(fr*o.gvc(X,1))).*abs(sin(fr*o.gvc(X,2)));
cof = mean(K(:));
% imagesc(K); title('variable coefficient plot');

varlap = @(u) - fftn(o.div( o.scalevec( o.grad(ifftn(u) ),K)));


b= varlap(fftn(u));

mvec = @(x)  o.C( varlap( o.S(x)) );
pvec = @(x) -cof*o.C(fftn(o.inv_laplacian(ifftn(o.S(x)))));

TOL = 1E-6;
MAXIT = 10*N;

fprintf('Next we report gmres output and relative error in the solution\n');

fprintf('No preconditioner\n');
tic; y1 = gmres(  mvec , b(:), [], TOL, MAXIT); toc
norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

fprintf('Constant coefficient laplace preconditioner\n'); 
tic;y0 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
norm(ifftn(o.S(y0))-u,inf)/norm(u,inf)


if 1
  L0=getOperator(@(x) cof*o.C(lap(o.S(x))),u(:));
  L =getOperator(@(x) o.C(varlap(o.S(x))),u(:));

  PL0 = mypinv(full(L0),1E-8);
  maxL=max(abs(L(:)));
  lI = abs(L)/maxL >1e-3;  
  fprintf('kept %2.2f  of the elements\n', full(sum(sum(lI(:))))/N^4);
  SL=full(lI.*L);
  PL = mypinv(SL,1E-8);
  pvec = @(x)PL*x;
  fprintf('Sparsified preconditioner with 1e-3 tolerance\n');
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

% Low rank off diagonal preconditioner. 
  M = N^2/2;
  D11 = L(1:M,1:M); D22=L(M+1:end, M+1:end);
  D=[[D11,0*D11];[0*D22, D22]];
  O= L-D;
  [U,S,V]=svd(full(O));
  r=1:200; Ur=U(:,r); Sr=S(r,r); Vr=V(:,r); 
  Or = Ur*Sr*Vr';

  PL = mypinv(D + Or ,1E-3);
  pvec = @(x)PL*x;
  fprintf('Single level HSS preconditioner with rank %d\n',r(end));
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

  SD11 = SL(1:M,1:M); SD22=SL(M+1:end, M+1:end);
  SD=[[SD11,0*D11];[0*D22, SD22]];
  PL = mypinv(SD + Ur*Sr*Vr',1E-8);
  pvec = @(x)PL*x;
  fprintf('Single level HSS preconditioner with sparsified diagonal\n');
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)



  m = 1:M/2;
  Da11=D11(m,m);Da22=D11(M/2+m,M/2+m);
  Db11=D22(m,m);Db22=D22(M/2+m,M/2+m);
  DD=0*D;
  DD(m,m)=Da11; m=M/2+m;DD(m,m)=Da22; m=M/2+m;DD(m,m)=Db11; m=M/2+m;DD(m,m)=Db22;
  OO = D-DD;
  [U,S,V]=svd(full(OO));
  Ur=U(:,r); Sr=S(r,r); Vr=V(:,r); 
  OOr = Ur*Sr*Vr';
  PL = mypinv(DD+Or+OOr,1E-8);
  pvec = @(x)PL*x;
  fprintf('Two level HSS preconditioner\n');
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

  PL = mypinv(full(DD),1E-8);
  pvec = @(x)PL*x;
  fprintf('Just diagonal of the two level HSS preconditioner\n');
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)


  PL = mypinv(L0+Or+OOr,1E-8);
  pvec = @(x)PL*x;
  fprintf('Two level HSS preconditioner with constant coefficent Laplacian in diagonal\n');
  tic;y1 = gmres(  mvec, b(:), [], TOL, MAXIT,pvec);  toc
  norm(ifftn(o.S(y1))-u,inf)/norm(u,inf)

end
