function stokes2(n)
if nargin<1, n=64; end
dim=2;
% set lattice dimensions
N = n*ones(dim,1);
 o = calc2(N);

% create points
X=o.regulargrid;

% create right hand side
f{1} = -2*sin(o.gvc(X,1)).*sin(o.gvc(X,2));
f{2} = cos(2*o.gvc(X,1)).^4;
f = o.crv(f);

p = o.inv_laplacian( -o.div(f) );
gradp = o.grad(p);
rhs = f - gradp;
u = o.inv_laplacian_vec(rhs);

subplot(3,1,1),surf(o.gvc(u,1)), axis off, shading interp, view(2), colormap bone;
subplot(3,1,2),surf(o.gvc(u,2)), axis off, shading interp, view(2), colormap bone;
subplot(3,1,3),surf(p), axis off, shading interp, view(2), colormap bone;
