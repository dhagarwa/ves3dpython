function poisson(n,dim)
if nargin<1, n=64; end
if nargin<2, dim=2; end;

% set lattice dimensions
N = n*ones(dim,1);
if dim==2,o = calc2(N);end
if dim==3,o = calc3(N);end

% create points
X=o.regulargrid;

% create right hand side
f = -2*sin(o.gvc(X,1));

u = o.inv_laplacian(f);
