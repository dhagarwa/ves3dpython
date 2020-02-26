function status= monitor(t,om,flag,mydata)

persistent cnt
matfile_stepsize = 100;
pngfile_stepsize = 2;

status = 0;
if strcmp(flag, 'setup_monitor'), cnt = 0; return; end
if strcmp(flag,'done'), return; end;

o=mydata{1};

omi=om;

om = om';
om = o.V(om(1,:));
u = -o.curl( o.inv_laplacian_vec(om));

infnorm_om  = norm(omi,inf);
twonorm_om  = norm(omi)/length(omi);
twonorm_u   = norm(o.C(u))/length(o.C(u));
CFL = 1/max(o.N*norm(o.C(u),inf));

str=sprintf('t=%.4e, Loo(om) = %.4e, L2(om)=%.2e, L2(u)=%.2e Loo(div(om))=%.2e CFL=%.2e\n',...
        t(1), infnorm_om, twonorm_om, twonorm_u, norm(o.C(o.div(om))), CFL);
fprintf(str);    

if ~isempty(mydata{2})
  omex = feval(mydata{2}, t(1));
  err = omex - om;
  infnorm_err = norm(o.C(err),inf)/norm(o.C(omex),inf);
  fprintf('||om_exact - om||_inf/ ||omega_exact||_inf : %e\n', infnorm_err);
end 

omn = o.inner(u,u);
%omn = o.gvc(om,3);

if 1

omn = interp3(omn,1,'cubic');

clf(gcf,'reset');
%h=vol3d('cdata', omn);  view(2), axis off; pause(0.02);
%f = omn(omn(:,:,o.N(3)/2));
[f1, f2,f3] = o.get_planes(omn);

%f1 =  omn(:,:,ceil(size(omn,3)/2));
%f2 = reshape(omn(:,ceil(size(omn,2)/2),:),[size(omn,1),size(omn,3)]);
%f3 = reshape(omn(ceil(size(omn,1)/2),:,:),[size(omn,2),size(omn,3)]);


ha = tight_subplot(2,3,[.001 .001], [.001 .001],[.001 .001]);

axes(ha(1));surf(f1); shading interp; axis equal, axis tight; axis off; colormap cool; view(2);
axes(ha(4));contour(f1),axis equal, axis off;  pause(0.02); colormap cool
axes(ha(2));surf(f2); shading interp; axis equal,axis tight; axis off; colormap cool; view(2);
axes(ha(5));contour(f2),axis equal, axis off;  pause(0.02); colormap cool
title(str);
axes(ha(3));surf(f3); shading interp; axis equal, axis tight; axis off; colormap cool; view(2);
axes(ha(6));contour(f3),axis equal, axis off;  pause(0.02); colormap cool
pause(0.2)
end

if 0
pngfile=sprintf('euler%05d_N%03d.png',cnt,o.N(1));
datafile=sprintf('euler%05d_N%03d.mat',cnt,o.N(1));

if mod(cnt, pngfile_stepsize), saveas(gcf,pngfile,'png'); end
if mod(cnt, matfile_stepsize*10), save(datafile, 'om','t', '-mat'); end
end


cnt = cnt + 1;
