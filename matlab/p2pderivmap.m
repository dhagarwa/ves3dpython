function diff = p2pderivmap(p_in, p_out)
%p1 nodes u-v coordinates according to p2 patch and its partial derivatives  
%diff = [du2/du1, du2/dv1, dv2/du1, dv2/dv1] 
    
    p1 = p_in.numPatch; p2 = p_out.numPatch;
    u = p_in.u; v = p_in.v;
    
    if p1==1 && p2==1 
        diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1 == 1 && p2==2
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==1 && p2==3
        diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==1 && p2==4
        diff =[1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==1 && p2==5
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = (-cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==1 && p2==6
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = (-cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];   
        
    elseif p1==2 && p2==1 
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==2 && p2==2
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==2 && p2==3
        diff = [1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==2 && p2==4
        diff =[1 + 0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==2 && p2==5
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = (cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==2 && p2==6
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);

        d3 = (cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];        
    
    elseif p1==3 && p2==1 
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==3 && p2==2
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==3 && p2==3
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==3 && p2==4
        diff =[1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==3 && p2==5
        d1 = -cos(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        
        d3 = -(sin(v)./cos(u).*(1./(tan(u).^2.*sin(v).^2+1)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*cos(v).*(cos(u).^2./(sin(u).^2.*sin(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==3 && p2==6
        d1 = cos(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        
        d3 = -(sin(v)./cos(u).*(1./(tan(u).^2.*sin(v).^2+1)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*cos(v).*(cos(u).^2./(sin(u).^2.*sin(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];        
        
    elseif p1==4 && p2==1 
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==4 && p2==2
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
    elseif p1==4 && p2==3
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
    elseif p1==4 && p2==4
        diff =[1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==4 && p2==5
        d1 = cos(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        
        
        d3 = (sin(v)./cos(u).*(1./(tan(u).^2.*sin(v).^2+1)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*cos(v).*(cos(u).^2./(sin(u).^2.*sin(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==4 && p2==6
        d1 = -cos(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*cos(p_in.v).^2).^(0.5);
        
        d3 = (sin(v)./cos(u).*(1./(tan(u).^2.*sin(v).^2+1)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*cos(v).*(cos(u).^2./(sin(u).^2.*sin(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*sin(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];        

    elseif p1==5 && p2==1 
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = -(cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
        
    elseif p1==5 && p2==2
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        
        d3 = (cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
    elseif p1==5 && p2==3
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = (1./sin(u).*(cos(v).^2./(cot(u).^2+cos(v).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(cos(u).*tan(v).*(sin(u).^2.*cos(v).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
    elseif p1==5 && p2==4
        d1 = -cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = -sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = -(1./sin(u).*(cos(v).^2./(cot(u).^2+cos(v).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (cos(u).*tan(v).*(sin(u).^2.*cos(v).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==5 && p2==5
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];

    elseif p1==5 && p2==6
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];
        
 
     elseif p1==6 && p2==1 
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = -(cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
        
    elseif p1==6 && p2==2
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        
        d3 = (cos(v)./cos(u).*(1./(tan(u).^2.*cos(v).^2+1)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(sin(u).*sin(v).*(cos(u).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
    elseif p1==6 && p2==3
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        
        d3 = -(1./sin(u).*(cos(v).^2./(cot(u).^2+cos(v).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = (cos(u).*tan(v).*(sin(u).^2.*cos(v).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];
    elseif p1==6 && p2==4
        d1 = cos(p_in.u).*sin(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
        d2 = sin(p_in.u).*cos(p_in.v)./(1-sin(p_in.u).^2.*sin(p_in.v).^2).^(0.5);
           
        d3 = (1./sin(u).*(cos(v).^2./(cot(u).^2+cos(v).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        d4 = -(cos(u).*tan(v).*(sin(u).^2.*cos(v).^2./(sin(u).^2.*cos(v).^2+cos(u).^2)).^(0.5))./(sin(u).^2.*cos(v).^2+cos(u).^2).^0.5;
        diff = [d1 d2 d3 d4];

    elseif p1==6 && p2==5
        diff = [0*p_in.u,0*p_in.u,0*p_in.u,0*p_in.u];

    elseif p1==6 && p2==6
        diff = [1+0*p_in.u,0*p_in.u,0*p_in.u,1+0*p_in.u];
        
        
    end
    

    diff(diff>10^10) = 0;
    diff(diff<-10^10) = 0;

end