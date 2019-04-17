function [res, dtA] = residual_2(meshdata, u, ufree, boundaries)

% boundaries = 0: freestream
% boundaries = 1: use boundaries

E2N = meshdata.E2N;
V = meshdata.V;
Area = meshdata.area;
I2E = meshdata.I2E;
B2E = meshdata.B2E;
In = meshdata.In;
Bn = meshdata.Bn;

nelem = length(u);
nInt = length(In);
nBound = length(Bn);
CFL = 0.9;

graduA_x = zeros(nelem,4);
graduA_y = zeros(nelem,4);

for i = 1:nInt
    eL = I2E(i,1);
    eR = I2E(i,3);
    normal = In(i,:);
    uL = u(eL,:);
    uR = u(eR,:);

    faceL = I2E(i,2);
    nodes = E2N(eL,:);
    nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);

    uavg = (uL+uR)./2;
    cont_x = uavg*normal(1)*deltal;
    cont_y = uavg*normal(2)*deltal;
    graduA_x(eL,:) = graduA_x(eL,:) + cont_x/Area(eL);
    graduA_y(eL,:) = graduA_y(eL,:) + cont_y/Area(eL);
    graduA_x(eR,:) = graduA_x(eR,:) - cont_x/Area(eR);
    graduA_y(eR,:) = graduA_y(eR,:) - cont_y/Area(eR);
end


for i = 1:nBound
    eL = B2E(i,1);
    normal = Bn(i,:);
%     ub = ufree;
    ub = u(eL,:)';

    faceL = B2E(i,2);
    nodes = E2N(eL,:);
    nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);

    cont_x = ub'*normal(1)*deltal;
    cont_y = ub'*normal(2)*deltal;
    graduA_x(eL,:) = graduA_x(eL,:) + cont_x/Area(eL);
    graduA_y(eL,:) = graduA_y(eL,:) + cont_y/Area(eL);   
end

gradu_x = graduA_x;
gradu_y = graduA_y;

res = zeros(nelem, 4);
sdl = zeros(nelem,1);

centroid = zeros(nelem,2);

for n = 1:nelem
    nodes = E2N(n,:);
    x_total = V(nodes(1),1) + V(nodes(2),1) + V(nodes(3),1);
    y_total = V(nodes(1),2) + V(nodes(2),2) + V(nodes(3),2);
    centroid(n,:) = [x_total, y_total]/3;
end

for i = 1:nInt
    eL = I2E(i,1);
    eR = I2E(i,3);
    normal = In(i,:);
    uL = u(eL,:)';
    uR = u(eR,:)';

    faceL = I2E(i,2);
    nodes = E2N(eL,:);
    nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);
    
    xavg = (V(nodes(2),1) + V(nodes(1),1))/2;
    yavg = (V(nodes(2),2) + V(nodes(1),2))/2;
    xvec_edge = [xavg, yavg];
    xdifL = xvec_edge - centroid(eL,:);
    xdifR = xvec_edge - centroid(eR,:);
    
    urecL = zeros(4,1);
    urecR = zeros(4,1);
    
    for j = 1:4
        gradL = [gradu_x(eL,j), gradu_y(eL,j)];
        urecL(j) = uL(j) + dot(gradL,xdifL);
        
        gradR = [gradu_x(eR,j), gradu_y(eR,j)];
        urecR(j) = uR(j) + dot(gradR,xdifR);
    end
    
    [Fhat, s] = roeflux(urecL, urecR, normal);
    
    res(eL,:) = res(eL,:) + Fhat'*deltal;
    res(eR,:) = res(eR,:) - Fhat'*deltal;
    sdl(eL) = sdl(eL) + s*deltal;
    sdl(eR) = sdl(eR) + s*deltal;
end

for i = 1:nBound
% for i = 1:15
    eL = B2E(i,1);
    normal = Bn(i,:);

    faceL = B2E(i,2);
    nodes = E2N(eL,:);
    nodes(faceL) = [];
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    deltal = sqrt(dx^2+dy^2);
    xavg = (V(nodes(2),1) + V(nodes(1),1))/2;
    yavg = (V(nodes(2),2) + V(nodes(1),2))/2;
    xvec_edge = [xavg, yavg];
    
    urecL = zeros(4,1);
    
    for j = 1:4
        xdifL = xvec_edge - centroid(eL,:);
        gradL = [gradu_x(eL,j), gradu_y(eL,j)];
        urecL(j) = u(eL,j) + dot(gradL,xdifL);
    end
    
    if boundaries == 0
        [Fhat, s] = roeflux(urecL,ufree,normal);
    elseif boundaries == 1
        if abs(normal(2)) > 0
            [Fhat, s] = wallflux(urecL, normal);
        elseif normal(1) < 0
            [Fhat, s] = influx(urecL, normal);
        else
            [Fhat, s] = outflux(urecL, normal);
        end
    end
    
    res(eL,:) = res(eL,:) + Fhat'*deltal;
    sdl(eL) = sdl(eL) + s*deltal;
end

dtA = 2*CFL./sdl;

end