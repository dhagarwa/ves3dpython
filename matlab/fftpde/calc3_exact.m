classdef calc3_exact 
properties
  x
  y 
  z
end

methods


%/******************************************************/
function o = calc3_exact
  syms xx  yy zz real
  o.x = xx;
  o.y = yy;
  o.z = zz;
end

%/******************************************************/
function gf = gradf(o,f)
% gradient of a scalar
gf = [f; f; f;];
gf(1) = diff(f,o.x);
gf(2) = diff(f,o.y);
gf(3) = diff(f,o.z);
end

%/******************************************************/
function gv = gradv(o,v)
g1 = o.gradf(v(1))';
g2 = o.gradf(v(2))';
g3 = o.gradf(v(3))';

gv = [g1; g2; g3];
end

%/******************************************************/
function f = divv(o,v)
gv = o.gradv(v);
f = gv(1,1) + gv(2,2) + gv(3,3);
end

%/******************************************************/
function w = curl(o, v)
gv = o.gradv(v);
gv_x = gv(1,:);
gv_y = gv(2,:);
gv_z = gv(3,:);

wx= gv_z(2) - gv_y(3);
wy= gv_x(3) - gv_z(1);
wz= gv_y(1) - gv_x(2);

w = [wx; wy; wz];
end

%/******************************************************/
function dT = divT(o,T)
fx = o.divv(T(1,:));
fy = o.divv(T(2,:));
fz = o.divv(T(3,:));
dT = [fx;fy;fz;];
end

%/******************************************************/
function T=outer(o, v, w)
v=v(:);
w=w(:);
T = [v*w(1), v*w(2), v*w(3)];
end

%/******************************************************/
function f=inner(o, v, w);
f = v(1)*w(1) + v(2)*w(2) + v(3)*w(3);
end


%/******************************************************/
function f=gaussian(o,X0,sigma,p)
d = [o.x;o.y;o.z] - X0(:); d=sqrt(d'*d);
arg = d/sigma;
if nargin>2, arg = arg.^p; else arg = arg.*arg; end
f = 1/sigma * exp(-arg);
end

%/******************************************************/
function f=trig(o,K)
arg = [o.x;o.y;o.z]' * K(:);
f = exp( 1i*arg);
end


%/******************************************************/
function test(o)
syms al sigma real
X0=[0;0;0]; %sigma=0.2;
K =[0;0;2];
ps = [o.gaussian(X0,sigma,4);0;0*real(o.trig(K))];
%ps = [real(o.trig(K));0;0]
u  = o.curl(ps);
om = o.curl(u);

syms tt real;
syms b real;
at = cos(b*pi*2*tt);
u = simple(at*u)
om = simple(o.curl(u))

Q= o.outer(om, u);
Q= Q-Q';
rhs = diff(om,tt) + o.divT(Q)

end

function bellmarkus(o)
  ro = 0.15; de =0.0333; ep=0.05; be = 15;
   x = o.x;
   y = o.y;
   z = o.z;
   U{1} = tanh( (ro - sqrt(z*z + y*y))/de);
   U{2} = 0*x;
   U{3} = ep * exp( -be*(x*x + y*y) );

   u = [U{1};U{2};U{3}]
   om = o.curl(u)
   o.divv(u)
   o.divv(om)
end


end % methods
end % class
