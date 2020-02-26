function status= monitor2(t,om,flag,mydata)

persistent cnt
matfile_stepsize = 100;
pngfile_stepsize = 2;

status = 0;
if strcmp(flag, 'setup_monitor'), cnt = 0; return; end
if strcmp(flag,'done'), return; end;
o=mydata{1};
omi=om;
om = om';
om = o.S(om(1,:));
u = -o.gradt( o.inv_laplacian(om));
infnorm_om  = norm(omi,inf);
twonorm_om  = norm(omi)/length(omi);
twonorm_u   = norm(o.C(u))/length(o.C(u));
CFL = 1/max(o.N*norm(o.C(u),inf));
str=sprintf('t=%.4e, Loo(om) = %.4e, L2(om)=%.2e, L2(u)=%.2e CFL=%.2e\n',...
        t(1), infnorm_om, twonorm_om, twonorm_u,CFL);
fprintf(str);    

if ~isempty(mydata{2})
  omex = feval(mydata{2}, t(1));
  err = omex - om;
  infnorm_err = norm(o.C(err),inf)/norm(o.C(omex),inf);
  fprintf('||om_exact - om||_inf/ ||omega_exact||_inf : %e\n', infnorm_err);
end 

if mod(cnt, pngfile_stepsize)
    surf(om),axis off, shading interp, view(2), colormap bone;
    pause(0.05);
end

omn = o.inner(u,u);

cnt = cnt + 1;
