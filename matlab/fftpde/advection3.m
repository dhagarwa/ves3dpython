clear all; clear globals
dim=3;

Nx = 96;
Ny = Nx;
Nz = Nx;
N = [Nx, Ny, Nz];
freq = [3,2,2];

o = calc3([Nx,Ny,Nz]);
o.use_real = true;

%om = calc3_functions(o, 'trigonometric', 'freq',freq);
om =  calc3_functions(o,'gaussian','sigma',pi/6); %Nx>=96 for good accuracy
u  = o.const_vec([1,1,1]);
rhs = @(om) -o.div( o.scalevec( u, om ) ) ;

% o.C() and o.V() are used to convert between vectors and tensor format used in calc3
[T,Omt]=ode45(@(t,om) o.C(rhs(o.S(om) )) ,[0,1],o.C(om));



