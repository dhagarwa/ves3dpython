clear all; clear globals;
dim=3;
Nx = 96; 
Ny=Nx; 
Nz=Nx;

freq = [3,2,2];
comp=2;

o = calc3([Nx, Ny, Nz]');
%o.use_hou_filtering;
o.use_twothirds_filtering;
% exact function  exp( i K.X)
X=o.regulargrid;


[fc,grad_fc] = calc3_functions(o, 'trigonometric', 'freq',freq);
%[fc,grad_fc] = calc3_functions(o,'gaussian','sigma',pi/6); %Nx>=96 for good accuracy

f = real(fc);
grad_f = real(grad_fc);

% fft differentiation
ap_grad_f = o.grad(f);

% error 
err_grad_f = norm(o.C(grad_f-ap_grad_f))/norm(o.C(grad_f))

str = 'either bug of aliasing errors';

% norm of error

assert( err_grad_f < 1e-13,str);
disp('gradient ok');

%% test vector divergence and curl

%create vector;
[vx, gvx]= calc3_functions(o,'trigonometric', 'freq',[2,1,3]);
[vy, gvy]= calc3_functions(o,'trigonometric', 'freq',[1,2,4]);
[vz, gvz]= calc3_functions(o,'trigonometric', 'freq',[4,3,1]);
tmp{1}=vx; tmp{2}=vy; tmp{3}=vz;
V = o.crv(tmp);

ap_divV = o.div(V);
divV= o.gvc(gvx,1) + o.gvc(gvy,2) + o.gvc(gvz,3);
assert(norm( o.C( ap_divV - divV)) <1e-10,str);
disp('div ok');

%%
ap_curlV = o.curl(V);
tmp{1}=o.gvc(gvz, 2) - o.gvc( gvy, 3);                                                        	
tmp{2}=o.gvc(gvx, 3) - o.gvc( gvz, 1);                                                          
tmp{3}=o.gvc(gvy, 1) - o.gvc( gvx, 2);
curlV = o.crv( tmp );     
assert( norm( o.C( ap_curlV - curlV)) <1e-10,str);
disp('curl ok');

%% check laplacian;
rho = o.div(gvy);
vy_ap = o.inv_laplacian(rho);
assert(norm(o.C(vy_ap-vy))<1e-10,str);
disp('inverse laplacian ok');

%% check inverse vector laplacian
rho = o.div_T(o.grad_V(V));
Vap = o.inv_laplacian_vec(rho);
assert(norm(o.C(Vap - V))<1e-10,str);
disp('inverse vec laplacian ok');



%%
%{ JUNK[
% function defined using symoblic toolbox
%syms xs ys zs;
%fs = imag(exp( 1i* (freq(1)*xs + freq(2)*ys + freq(3)*zs) ));
%fs_v = subs(fs, {xs,ys,zs}, {X(:,:,:,1),X(:,:,:,2),X(:,:,:,3)});
%grad_fs = [diff(fs,xs); diff(fs,ys); diff(fs,zs);];
%grad_fs_v = subs(grad_fs(comp), {xs,ys,zs}, {X(:,:,:,1),X(:,:,:,2),X(:,:,:,3)});

% errors between analytic
%col=@(x)x(:);
%norm(col(f-fs_v))
%norm(col(grad_fs_v-grad_f(:,:,:,comp) ))
%subplot(2,1,1), surf(real(reshape(ap_grad_f(:,:,2,2),[o.N(2),o.N(1)]))), shading interp; 
%subplot(2,1,2), surf(real(reshape(grad_f(:,:,2,2),[o.N(2),o.N(1)]))), shading interp; 


%}




