function [] = testPoissonSolverMat()
    clc
    

    a = 0;
    b = 2*pi;
    N = 64; % number of cells; N+1 is the number of grid points
    L = b-a; % size of the domain
    h = L/N; % grid spacing
    x = linspace(a,b,N+1)'; % column vector of grid points

    phi_f = @(t) t.*cos(t); % inline function for the exact solution
    rho_f = @(t) 2.*sin(t) + t.*cos(t); % inline function for the exact right-hand-side

    phi_exact = phi_f(x); % exact solution at the grid points

    alpha = phi_f(a); % boundary condition
    beta = phi_f(b); % boundary condition

    phi = zeros(N-1,1); % column vector for the solution
    A = zeros(N-1,N-1); % matrix of discrete Laplacian.

    rho = rho_f(x(2:N)); % exact rho at the grid points
    rho(1) = rho(1) + alpha/h^2; % add boundary term
    rho(N-1) = rho(N-1) + beta/h^2; % add boundary term

    dA = diag( 2*ones(1,N-1) ); % diagonal matrix
    dAp1 = diag( -1*ones(1,N-2), 1 ); % super-diagonal matrix
    dAm1 = diag( -1*ones(1,N-2), -1 ); % sub-diagonal matrix
    A = (dA + dAp1 + dAm1);
    
    A = A/h^2;

    eigv = eig(A) % for curiosity, we check the eigenvalues of A
    cond_num = cond(A)
    
    S = svd(A);
    plot(1:size(S,1),S);
    
    
    phi = A\rho; % solving the linear system

    solver_err = max( abs(A*phi - rho) ) % Did the solver do his job?
end