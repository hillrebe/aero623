function [Rvec, dtA] = residual_1(meshdata, u, ufree, boundaries)

E2N = meshdata.E2N;
V = meshdata.V;
% Area = meshdata.area;
I2E = meshdata.I2E;
B2E = meshdata.B2E;
In = meshdata.In;
Bn = meshdata.Bn;

nelem = length(u);
nInt = length(In);
nBound = length(Bn);
CFL = 0.95;

sdl = zeros(nelem,1);
Rvec = zeros(nelem,4);

localnodes = [1 2 3 4 1];

for i = 1:nInt
    eL = I2E(i,1);
    eR = I2E(i,3);
    normal = In(i,:);
    uL = u(eL,:)';
    uR = u(eR,:)';
    [Fhat, s] = roeflux(uL,uR,normal);

    faceL = I2E(i,2);
%     faceR = I2E(i,4);
    nodes_k = localnodes(faceL:faceL+1);
    nodes = E2N(eL,nodes_k);
%     nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);
    sdl(eL) = sdl(eL) + s*deltal;
    sdl(eR) = sdl(eR) + s*deltal;

    Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal;
    Rvec(eR,:) = Rvec(eR,:) - Fhat'*deltal;
end

for i = 1:nBound
    eL = B2E(i,1);
    normal = Bn(i,:);
    uL = u(eL,:)';
    bgroup = B2E(i,3);

    if boundaries == 0
        [Fhat, s] = roeflux(uL,ufree,normal);
    elseif boundaries == 1
        if bgroup == 1
            [Fhat, s] = wallflux(uL, normal);
        else
            [Fhat, s] = roeflux(uL, ufree, normal);
        end
    end           

    faceL = B2E(i,2);
    nodes_k = localnodes(faceL:faceL+1);
    nodes = E2N(eL,nodes_k);
%     nodes = E2N(eL,:);
%     nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);
    sdl(eL) = sdl(eL) + s*deltal;

    Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal;
end

dtA = 2*CFL./sdl;
        
end