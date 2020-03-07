function val = fftinterp1(x, f)
%Function to interpolate using ffts at non grid points given in x
%1D function
%f is the function values on [-pi, pi] uniform periodic grid
    
   n = size(f, 1);
   x0 = (0:n-1)* 2*pi/n; 
   x0 = -pi + x0';
   %f = ftest(x0);
   f_hat =  fft(f); % TODO pass the fhats instead to avoid this calculation everytime
   k = [0:n/2 -n/2+1:-1]; k = k';
   %f_app = ifft(f_hat);   
%    x = (0:4*n-1)*2*pi/(4*n);
%    x = -pi + x';
%    f = ftest(x);
   
   val = zeros(size(x, 1), 1);
   for ii=1:size(x, 1)
       %size(exp(1i*k*x(ii)).*f_hat)
      val(ii) = sum(exp(1i*k*(x(ii)+pi)).*f_hat/n);  % matlab assumes 0, 2*pi interval
   end
%    plot(x, f, 'o');
%    hold on;
%    plot(x, f_app, 'x');
%    norm(f-f_app)

end

% function val = ftest(x)
%     val = exp(sin(x).^2 + cos(exp(cos(x).^5)) ) ;
% 
% end