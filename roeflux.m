function [Fhat, s] = roeflux(uleft, uright, n)

gamma = 1.4;
uvec = [uleft, uright];
rhoL = uvec(1,1);
rhoR = uvec(1,2);

d_rho = rhoR-rhoL;
d_rhou = uvec(2,2)-uvec(2,1);
d_rhov = uvec(3,2)-uvec(3,1);
d_rhoE = uvec(4,2)-uvec(4,1);

for i = 1:2
    rho = uvec(1,i);
    rhou = uvec(2,i);
    rhov = uvec(3,i);
    rhoE = uvec(4,i);
    
    utemp = rhou/rho;
    vtemp = rhov/rho;
    q2 = utemp^2 + vtemp^2;
    p = (gamma-1)*(rhoE-0.5*rho*q2);
    Etemp = rhoE/rho;
    H = Etemp + p/rho;
    c = sqrt(gamma*p/rho);
    
    Fvec_x(1:4,i) = [rhou; rhou*utemp + p; rhou*vtemp; rhou*H];
    Fvec_y(1:4,i) = [rhov; rhov*utemp; rhov*vtemp + p; rhov*H];
    
    vars(1:8,i) = [uvec(:,i); utemp; vtemp; H; c];
end

uRoe = roeavg(rhoL,rhoR, vars(5,1),vars(5,2));
vRoe = roeavg(rhoL,rhoR, vars(6,1),vars(6,2));
HRoe = roeavg(rhoL,rhoR, vars(7,1),vars(7,2));
% cRoe = roeavg(rhoL,rhoR, vars(8,1),vars(8,2));

vvecRoe = [uRoe vRoe];
qRoe = norm(vvecRoe);
cRoe = sqrt((gamma-1)*(HRoe-0.5*qRoe^2));

udot = dot(vvecRoe,n);

eps = 0.1*cRoe;
lambda = [udot+cRoe; udot-cRoe; udot];
for i = 1:3
    if abs(lambda(i)) < eps
        lambda(i) = (eps^2+lambda(i)^2)/(2*eps);
    end
end

lamb1 = lambda(1);
lamb2 = lambda(2);
lamb3 = lambda(3);

s1 = 0.5*(abs(lamb1)+abs(lamb2));
s2 = 0.5*(abs(lamb1)-abs(lamb2));

G1 = (gamma-1)*(qRoe^2/2*d_rho - dot(vvecRoe,[d_rhou d_rhov]) + d_rhoE);
G2 = -udot*d_rho + dot([d_rhou d_rhov],n);

C1 = (G1/cRoe^2)*(s1-abs(lamb3)) + G2*s2/cRoe;
C2 = G1*s2/cRoe + (s1-abs(lamb3))*G2;

for i = 1:4
    for j = 1:2
        F_ind = [Fvec_x(i,j) Fvec_y(i,j)];
        Ftemp(i,j) = dot(F_ind,n);
    end
end

term1 = 0.5*(Ftemp(:,1)+Ftemp(:,2));
term2 = 0.5 * [abs(lamb3)*d_rho + C1;
    abs(lamb3)*d_rhou + C1*vvecRoe(1) + C2*n(1);
    abs(lamb3)*d_rhov + C1*vvecRoe(2) + C2*n(2);
    abs(lamb3)*d_rhoE + C1*HRoe + C2*udot];
Fhat = term1-term2;

s = abs(udot) + cRoe;

end