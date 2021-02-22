function sdiff = p2psecondderivmap(p_in, p_out)
%p1 nodes u-v coordinates according to p2 patch and its partial derivatives  
%first_diff = [du2/du1, du2/dv1, dv2/du1, dv2/dv1] 
%diff = [d^2u2/du1^2 d^2u2/du1dv1, d^2u2/du1dv1 d^2u2/d^2v1, d^2v2/d^2u1 d^2v2/du1dv1, d^2v2/dv1du1 d^2v2/dv1^2  ]
    
    p1 = p_in.numPatch;
    p2 = p_out.numPatch;
    u1 = p_in.u + 10^8*eps; v1 = p_in.v + 10^8*eps;
    syms u v;
    
    if p1==1 && p2==1 
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
        %sdiff = [0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u];
    elseif p1 == 1 && p2==2
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
        %sdiff = [0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u];
    elseif p1==1 && p2==3
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
        %sdiff = [0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u];
    elseif p1==1 && p2==4
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff =[1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
        %sdiff = [0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u, 0*p_in.u];
    elseif p1==1 && p2==5
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2 )^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2 )^(0.5);
        
        d3 = (-cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
        %sd1 = (cos(u).^2.*sin(u).*sin(v).^3)/(1 - sin(u).^2.*sin(v).^2).^(3/2) - (sin(u).*sin(v))./(1 - sin(u).^2.*sin(v).^2).^(1/2);
        %sd2 = (cos(u).*cos(v))./(1 - sin(u).^2.*sin(v).^2).^(1/2) + (cos(u).*cos(v).*sin(u).^2.*sin(v).^2)/(1 - sin(u).^2.*sin(v).^2).^(3/2);
        %sd3 = (cos(u)*cos(v))/(1 - sin(u)^2*sin(v)^2)^(1/2) + (cos(u)*cos(v)*sin(u)^2*sin(v)^2)/(1 - sin(u)^2*sin(v)^2)^(3/2)
        
    elseif p1==1 && p2==6
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = (-cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];   
        
    elseif p1==2 && p2==1 
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==2 && p2==2
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==2 && p2==3
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==2 && p2==4
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff =[1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==2 && p2==5
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = (cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];

    elseif p1==2 && p2==6
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);

        d3 = (cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];        
    
    elseif p1==3 && p2==1 
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==3 && p2==2
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==3 && p2==3
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==3 && p2==4
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff =[1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==3 && p2==5
        d1 = -cos(u)*cos(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        d2 = sin(u)*sin(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        
        d3 = -(sin(v)/cos(u)*(1/(tan(u)^2*sin(v)^2+1))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*cos(v)*(cos(u)^2/(sin(u)^2*sin(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];

    elseif p1==3 && p2==6
        d1 = cos(u)*cos(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        d2 = -sin(u)*sin(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        
        d3 = -(sin(v)/cos(u)*(1/(tan(u)^2*sin(v)^2+1))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*cos(v)*(cos(u)^2/(sin(u)^2*sin(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];        
        
    elseif p1==4 && p2==1 
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==4 && p2==2
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==4 && p2==3
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==4 && p2==4
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff =[1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==4 && p2==5
        d1 = cos(u)*cos(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        d2 = -sin(u)*sin(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        
        
        d3 = (sin(v)/cos(u)*(1/(tan(u)^2*sin(v)^2+1))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*cos(v)*(cos(u)^2/(sin(u)^2*sin(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];

    elseif p1==4 && p2==6
        d1 = -cos(u)*cos(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        d2 = sin(u)*sin(v)/(1-sin(u)^2*cos(v)^2)^(0.5);
        
        d3 = (sin(v)/cos(u)*(1/(tan(u)^2*sin(v)^2+1))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*cos(v)*(cos(u)^2/(sin(u)^2*sin(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*sin(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];        

    elseif p1==5 && p2==1 
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = -(cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
        
    elseif p1==5 && p2==2
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        
        d3 = (cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
    elseif p1==5 && p2==3
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = (1/sin(u)*(cos(v)^2/(cot(u)^2+cos(v)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(cos(u)*tan(v)*(sin(u)^2*cos(v)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
    elseif p1==5 && p2==4
        d1 = -cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = -sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = -(1/sin(u)*(cos(v)^2/(cot(u)^2+cos(v)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (cos(u)*tan(v)*(sin(u)^2*cos(v)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];

    elseif p1==5 && p2==5
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==5 && p2==6
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
        
 
     elseif p1==6 && p2==1 
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = -(cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
        
    elseif p1==6 && p2==2
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        
        d3 = (cos(v)/cos(u)*(1/(tan(u)^2*cos(v)^2+1))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(sin(u)*sin(v)*(cos(u)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
    elseif p1==6 && p2==3
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        
        d3 = -(1/sin(u)*(cos(v)^2/(cot(u)^2+cos(v)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = (cos(u)*tan(v)*(sin(u)^2*cos(v)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];
    elseif p1==6 && p2==4
        d1 = cos(u)*sin(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
        d2 = sin(u)*cos(v)/(1-sin(u)^2*sin(v)^2)^(0.5);
           
        d3 = (1/sin(u)*(cos(v)^2/(cot(u)^2+cos(v)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        d4 = -(cos(u)*tan(v)*(sin(u)^2*cos(v)^2/(sin(u)^2*cos(v)^2+cos(u)^2))^(0.5))/(sin(u)^2*cos(v)^2+cos(u)^2)^0.5;
        %diff = [d1 d2 d3 d4];

    elseif p1==6 && p2==5
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
        %diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];

    elseif p1==6 && p2==6
        d1 = 1; d2 = 0; d3 = 0; d4 = 1;
        %diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
        
        
    end
    
    
%diff = [d^2u2/du1^2 d^2u2/du1dv1, d^2u2/du1dv1 d^2u2/d^2v1, d^2v2/d^2u1 d^2v2/du1dv1, d^2v2/dv1du1 d^2v2/dv1^2  ]    
    sd1 = diff(d1, u);
    sd2 = diff(d1, v);
    sd3 = diff(d2, u);
    sd4 = diff(d2, v);
    sd5 = diff(d3, u);
    sd6 = diff(d3, v);
    sd7 = diff(d4, u);
    sd8 = diff(d4, v);
    
    
    sd1_ = double(subs(sd1, {'u' 'v'}, {u1 v1}));
    sd2_ = double(subs(sd2, {'u' 'v'}, {u1 v1}));
    sd3_ = double(subs(sd3, {'u' 'v'}, {u1 v1}));
    sd4_ = double(subs(sd4, {'u' 'v'}, {u1 v1}));
    sd5_ = double(subs(sd5, {'u' 'v'}, {u1 v1}));
    sd6_ = double(subs(sd6, {'u' 'v'}, {u1 v1}));
    sd7_ = double(subs(sd7, {'u' 'v'}, {u1 v1}));
    sd8_ = double(subs(sd8, {'u' 'v'}, {u1 v1}));

    sdiff = [sd1_, sd2_, sd3_, sd4_, sd5_, sd6_, sd7_, sd8_];
    sdiff(sdiff>10^10) = 0;
    sdiff(sdiff<-10^10) = 0;

end