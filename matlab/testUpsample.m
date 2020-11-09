function [] = testUpsample()
    %testing the interpolation - upsample + filtering 
    m = 15; n =15;
    [u,v]=ndgrid(2*pi*(1:m)/(m+1),  (2*pi)*(1:n)/(n+1));
    
    
    %f = sin(u).*cos(v);
    %stem3(f);
    u = u(:);
    v = v(:);
    f = sin(7*u).*sin(v); % + sin(10*u).*sin(v) + sin(5*u).*sin(6*v);
    f_int = upsample2(f, []); %interpolated f
    [u_up,v_up]=ndgrid(2*pi*(0:2*m+1)/(2*m+2),  (2*pi)*(0:2*n+1)/(2*n+2));
    u_up = u_up(:);
    v_up = v_up(:);
    f_up = sin(7*u_up).*sin(v_up);% + sin(10*u_up).*sin(v_up) + sin(5*u_up).*sin(6*v_up);    
    
    error = max(abs(f_int - f_up))
    
    
end