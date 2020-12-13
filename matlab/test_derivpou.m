function [] = test_derivpou()
%test resolution of derivative of pou function 
%see how any frequncies resolve it accurately


      clc

    m = 127;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S1 = Surface(patches, [0, 0, 0]);
    p1 = patches(1);
    %pou_err = checkerror(sin(16*p.u).*sin(16*p.v),p)
    p1.pou = pou(p1, p1.numPatch, p1.r);
    [pou1_du, pou1_dv] = deriv_pou(p1, p1.numPatch, p1.r);
%     figure;
%     stem3(reshape(p1.pou,[p1.Nu,p1.Nv]));
     figure;
     stem3(reshape(pou1_du,[p1.Nu,p1.Nv]));
%     figure;
%     stem3(reshape(pou1_dv,[p1.Nu,p1.Nv]));
    
    [pou1_du_app, pou1_dv_app] = p1.grad_patch(p1.pou);
    
%     figure;
%     stem3(reshape(pou1_du_app,[p1.Nu,p1.Nv]));
%     figure;
%     stem3(reshape(pou1_dv_app,[p1.Nu,p1.Nv]));
%     
    
    err_du = max(abs(pou1_du_app-pou1_du))
    err_dv = max(abs(pou1_dv_app-pou1_dv))
    
    m = 2*(m+1)-1;
    n = m;
    R = 1;
    
    %theta = pi/2;
    %phi = pi/3;
    %trg = [R*sin(theta)*cos(phi) R*sin(theta)*sin(phi) R*cos(theta)]
    trg = [0.997425 -0.0507148 -0.0507148];
    norm_trg = sqrt(norm(trg));
    patches = [];
    for i=1:6
       patch =  standardSpherePatch(m, n, i, R);
       
       patch.q_sl = fooVec(patch.r);
       patch.q_dl = fooVec(patch.r);
       patches = [patches patch];
       %val = val + DLSmoothPatch(trg, patch, f);
        
    end
    S2 = Surface(patches, [0, 0, 0]);
    p2 = patches(1);
    p2.pou = pou(p2, p2.numPatch, p2.r);
    
  
    
    pou_up = upsample2(p1.pou,p1);
    %figure;
    %stem3(reshape(pou_up,[p2.Nu,p2.Nv]));
    
    
    
    
    err_pou = max(abs(pou_up-p2.pou))
    %max(abs(p1_pou_up-p2.pou))
    
    






end


function val = pou(obj, patch, r) %patch is patch num here
    [v1, u1, rho1] = cart2sph(r(:,1), r(:,2), r(:,3));
    %u = pi/2 - u;
    %u; v;
    if patch==1
        r0 = [0,obj.R,0];
    elseif patch ==2    
        r0 = [0,-obj.R,0];
    elseif patch==5
        r0 = [0,0,obj.R];
    elseif patch==6
        r0 = [0,0,-obj.R];
    elseif patch==3
        r0 = [obj.R,0,0];
    elseif patch==4
        r0 = [-obj.R,0,0];
        
    end
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2-u0;
    %u0; v0;
    d = (5/12) * pi*obj.R;
    nn = size(r, 1);
    val = zeros(nn, 1);
    t_vec = zeros(nn, 1);
    for ii=1:nn
        t = greatCircleDistance(u0, v0, u1(ii), v1(ii), obj.R)/d;
        t_vec(ii) = t;
        if t>=1
            val(ii) = 0;
        elseif t ==0
            val(ii) = 1;
        else
            val(ii) = exp((2*exp(-1/t))/ (t-1)); %pou from bruno paper
            %val(ii) = exp((2*exp(-1/t))/ (t-1)); %derivative of pou from bruno paper
            %val(ii) = exp((1*exp(-1/t))/ (t-1)); %pou from bruno paper
            %val(ii) = exp((2*exp(-1/t^2))/ (t^2-1)); %MORE SMOOTH pou 
            %val(ii) = exp(t^2/ (t^2-1)); %POU from beale paper
        end
    end
    
    %plot(val)
    %disp("plotting great circle distance");
end

function [du, dv] = deriv_pou(obj, patch, r) %patch is patch num here
    [v1, u1, rho1] = cart2sph(r(:,1), r(:,2), r(:,3));
    %u = pi/2 - u;
    %u; v;
    if patch==1
        r0 = [0,obj.R,0];
    elseif patch ==2    
        r0 = [0,-obj.R,0];
    elseif patch==5
        r0 = [0,0,obj.R];
    elseif patch==6
        r0 = [0,0,-obj.R];
    elseif patch==3
        r0 = [obj.R,0,0];
    elseif patch==4
        r0 = [-obj.R,0,0];
        
    end
    [v0, u0, rho0] = cart2sph(r0(1), r0(2), r0(3));
    %u0 = pi/2-u0;
    %u0; v0;
    d = (5/12) * pi*obj.R;
    nn = size(r, 1);
    du = zeros(nn, 1);
    dv = zeros(nn, 1);
    t_vec = zeros(nn, 1);
    for ii=1:nn
        t = greatCircleDistance(u0, v0, u1(ii), v1(ii), obj.R)/d;
        t_vec(ii) = t;
        if t>=1
            du(ii) = 0;
            dv(ii) = 0;
        elseif t ==0
            du(ii) = 0;
            dv(ii) = 0;
        else
            %val(ii) = exp((2*exp(-1/t))/ (t-1)); %pou from bruno paper
            phi = exp((2*exp(-1/t))/ (t-1)); %derivative of pou from bruno paper
            dt = 2*(exp(-1/t)/((t-1)*t^2) - exp(-1/t)/(t-1)^2);
            dphi_dt = phi*dt;
            u = obj.u(ii); v = obj.v(ii);
            dt_du = -cos(u)*sin(v)/(d*(1-sin(u)^2*sin(v)^2)^(0.5));
            dt_dv = -cos(v)*sin(u)/(d*(1-sin(u)^2*sin(v)^2)^(0.5));
            du(ii) = dphi_dt*dt_du;
            dv(ii) = dphi_dt*dt_dv;
            %val(ii) = exp((1*exp(-1/t))/ (t-1)); %pou from bruno paper
            %val(ii) = exp((2*exp(-1/t^2))/ (t^2-1)); %MORE SMOOTH pou 
            %val(ii) = exp(t^2/ (t^2-1)); %POU from beale paper
        end
    end
    
    %plot(val)
    %disp("plotting great circle distance");
end