A = rand(30000, 30000);
B = gpuArray(A);
tic 
Ahat = fft(A);
t1 = toc
tic;
Bhat = fft(B);
t2 = toc
speedup = t1/t2 
