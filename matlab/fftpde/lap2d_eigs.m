clear all; clear globals;

N=16;  % mesh size

% create Fourier object
o=calc2([N,N]);   

% consp
id=eye(N*N);
A=zeros(N^2);
for i=1:N^2
  x=o.laplacian(reshape(id(i,:),N,N));
  A(:,i)=-x(:);
end

% symmetrize to avoid eigenvalue errors
A=1/2 * (A+A');

[Q,L]=eig(A);

any(imag(diag(L)))
any((-diag(L))>0)
semilogy(-diag(L))
max(abs(imag(diag(L)))
max(abs(imag(diag(L))
max(abs(imag(diag(L))))
norm(A-A')/norm(A)
max(abs(imag(diag(L))))
max(abs(real(diag(L))))
for i=1:N^2; surface(reshape(Q(:,ind(i)),N,N)); view(2); shading interp; pause; clf;end
[~,ind]=sort(diag(L),'ascend');
[~,ind]=sort(real(diag(L)),'ascend');
for i=1:N^2; surface(reshape(Q(:,ind(i)),N,N)); view(2); shading interp; pause; clf;end
for i=1:N^2; surface(reshape(real(Q(:,ind(i)),N,N))); view(2); shading interp; pause; clf;end
for i=1:N^2; surface(reshape(real(Q(:,ind(i))),N,N)); view(2); shading interp; pause; clf;end
[~,ind]=sort(diag(L),'ascend');
[~,ind]=sort(real(abs(diag(L))),'ascend');
for i=1:N^2; surface(reshape(real(Q(:,ind(i))),N,N)); view(2); shading interp; pause; clf;end
ind
