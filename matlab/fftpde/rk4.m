function y=rk4(dydt,tint,y0,ops)

dt = ops.InitialStep;
T0=tint(1);
T1=tint(2);
Tc=T0;
M=(T1-T0)/dt;  
M=round(M);

mon = ops.OutputFcn;
mon(tint(1),y0,'setup_monitor');
mon(tint(1),y0,[]);

for j=1:M
  F0 = dydt(Tc,y0);
  F1 = dydt(Tc+dt/2, y0+ dt/2*F0);
  F2 = dydt(Tc+dt/2, y0+ dt/2*F1);
  F3 = dydt(Tc+dt, y0+ dt*F2);
  y= y0 + dt/6*(F0+2*F1+2*F2+F3);

  mon(Tc+dt , y, []);
  y0=y;
  Tc=Tc+dt;
end

mon(tint(1)+M*dt,y0,'done');
