function [cl, cpsort] = calcOutputs(umat, meshdata, plotswitch)

% plotswitch = 1 to plot, 0 to not plot

Area = meshdata.Area;
B2E = meshdata.B2E;
Bn = meshdata.Bn;
V = meshdata.V;
E2N = meshdata.E2N;

pinf = 1;
Minf = 0.2;
h = 1;
gamma = 1.4;
Tt = 1+(gamma-1)/2*Minf^2;
pt = Tt^(gamma/(gamma-1));
R = 1;
rhot = pt/(R*Tt);
st = pt/(rhot^gamma);

% find edges along bottom boundary
B2Eairfoil = B2E(B2E(:,3)==1,:);
Bnairfoil = Bn(B2E(:,3)==1,:);
nAirfoil = length(Bnairfoil);

% drag and lift coefficients
cd_num = 0;
cl_num = 0;
% ent_num = 0;
% ar = 0;
localnodes = [1 2 3 4 1];

for i = 1:nAirfoil
    elem = B2Eairfoil(i,1);
    elemvec(i) = elem;
    face = B2Eairfoil(i,2);
    nodes_k = localnodes(face:face+1);
    nodes = E2N(elem,nodes_k);
    dx = V(nodes(2),1) - V(nodes(1),1);
    dy = V(nodes(2),2) - V(nodes(1),2);
    dl = sqrt(dx^2+dy^2);
    normal = Bnairfoil(i,:);
    
    rho = umat(elem,1);
    rhou = umat(elem,2);
    rhov = umat(elem,3);
    rhoE = umat(elem,4);

    utemp = rhou/rho;
    vtemp = rhov/rho;
    q2 = utemp^2 + vtemp^2;
    p = (gamma-1)*(rhoE-0.5*rho*q2);
    
    cl_num = cl_num + (p-pinf)*normal(2)*dl;
    cd_num = cd_num + (p-pinf)*normal(1)*dl;
    
    xmid(i) = (V(nodes(2),1) + V(nodes(1),1))/2;
%     x1 = V(nodes(2),1);
%     x2 = V(nodes(1),1);
%     xleft(i) = min([x1,x2]);
%     xright(i) = max([x1,x2]);
    cpnum = p-pinf;
    cpden = (gamma/2)*pinf*Minf^2;
    cp(i) = cpnum/cpden;
end

% xtemp = sort(xmid);
elemtemp = sort(elemvec);
for i = 1:nAirfoil
%     k = find(xmid==xtemp(i));
    k = find(elemvec == elemtemp(i));
    xsort(i) = xmid(k);
    cpsort(i,1) = cp(k);
end

% for i = 1:nLower
%     k = (i-1)*2 + [1:2];
%     xsort(k) = [xleft_temp(i); xright_temp(i)];
%     cpsort(k) = [cpsort_temp(i); cpsort_temp(i)];
% end

% nelem = length(E2N);

% for j = 1:nelem
%     elem = j;
%     rho = umat(elem,1);
%     rhou = umat(elem,2);
%     rhov = umat(elem,3);
%     rhoE = umat(elem,4);
% 
%     utemp = rhou/rho;
%     vtemp = rhov/rho;
%     q2 = utemp^2 + vtemp^2;
%     p = (gamma-1)*(rhoE-0.5*rho*q2);
%     s = p/rho^gamma;
%     sterm = ((s/st)-1)^2;
%     ent_num = ent_num + sterm*Area(j);
%     ar = ar + Area(j);
% end

den = (gamma/2)*pinf*Minf^2*h;
cl = cl_num/den;
cd = cd_num/den;
% Es = sqrt(ent_num/ar);
% cp = 1;

if plotswitch == 1
    figure
    plot(xsort,cpsort)
    set(gca,'Ydir','reverse')
    xlabel('x-position along bottom boundary')
    ylabel('c_p')
end

end