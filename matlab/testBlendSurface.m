function [] = testBlendSurface()

    %membrane properties
    kb = 0;
    Es = 2;
    Ed = 20;
    vesprops = [0;Es;Ed]; %kb, Es, Ed
    
    %flow properties
    shear_rate = 5;
    
    %Initialize a surface
    
    clc

    m = 15;
    n = 15;
    R = 1;
    

    %trg = [0.997425 -0.0507148 -0.0507148];
    %norm_trg = sqrt(norm(trg))
    patches = [];
    numPatches = 6;
    for i=1:numPatches
       patch =  standardSpherePatch(m, n, i, R);
       patches = [patches patch];
        
    end
    S = Surface(patches, vesprops);
    S.updateSurface(2*ones(size(S.getPosition(), 1), 3), 1); %move to (2, 2, 2)
    
    %makePosition(S, 1);
    S.updateStale();
    %plotSurface(S);
    r_before_blend = S.getPosition();
    S.blendSurface();
    r_after_blend = S.getPosition();
    size(r_after_blend)
    error_blend = norm(r_before_blend - r_after_blend)/(m*n*6)
    plot(1:225, r_before_blend(1:225, 1));
    hold on;
    %plot(1:1350, r_after_blend(:, 1));
    %figure;
    %plotSurface(S);


end

function [] = makePosition(S, val)
    for ii=1:S.numPatches
        patch = S.patches(ii);
        patch.r = ones(patch.numNodes, 3);
        patch.rb = patch.r;
    end

end