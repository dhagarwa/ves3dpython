% Test of the implementation of the 2D FDM which can be used for solving
% PDE
% C. Weng
% DLR, Berlin
% 1st version: 24-May-2017

% the function to be tested and the analytic solution to its derivative
fun = @(x,y) x.^3.*cos(pi*y);
dfundx = @(x,y) 3*x.^2.*cos(pi*y);
dfundy = @(x,y) x.^3.*-1*pi.*sin(pi*y);


% parameters
npx = 16;
npy = 16;
xVec = linspace(-1,1,npx);
yVec = linspace(-1,1,npy);
dx = diff(xVec([1 2]));
dy = diff(yVec([1 2]));
n = 1;  % derivative order
ooa = 8; % order of accuracy of the FDM
tic
[Dx, Dy] = getNonCompactFDmatrix2D(npx,npy,dx,dy,n,ooa);
toc



% generate function vector
[XX,YY] = meshgrid(xVec,yVec);
x = XX(:);
y = YY(:);
funVec = fun(x,y);
dfundxAna = dfundx(x,y);
dfundyAna = dfundy(x,y);

% apply the diff. matrix
dfundxNum = Dx*funVec;
dfundyNum = Dy*funVec;


% Abs. error
dfundxErr = abs(dfundxNum-dfundxAna);
dfundyErr = abs(dfundyNum-dfundyAna);

%% plot
%*****  plot the function
figure(1)
clf
surf(XX,YY,reshape(funVec,npy,npx))
xlabel('x'),ylabel('y'),zlabel(func2str(fun))

% dF/dx
figure(2)
subplot(211)
surf(XX,YY,reshape(dfundxErr,npy,npx))
xlabel('x'),ylabel('y'),zlabel('Abs. Error(dF/dx)')

% dF/dy
subplot(212)
surf(XX,YY,reshape(dfundyErr,npy,npx))
xlabel('x'),ylabel('y'),zlabel('Abs. Error(dF/dy)')













