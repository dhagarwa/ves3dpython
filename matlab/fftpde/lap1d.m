%clear all; clear globals;

N=16;  % mesh size

% create Fourier object
o=calc1([N]);   

% consp
id=eye(N);
A=zeros(N);
for i=1:N
  x=o.laplacian(id(:,i));
  A(:,i)=-x(:);
end

% symmetrize to avoid eigenvalue errors
A=1/2 * (A+A');

[Q,L]=eig(A);

