function [Fb, s] = wallflux(uint, n)

n = n';
gamma = 1.4;
rhoplus = uint(1);
rhou = uint(2);
rhov = uint(3);
rhoE = uint(4);
utemp = rhou/rhoplus;
vtemp = rhov/rhoplus;
q2 = utemp^2 + vtemp^2;
p = (gamma-1)*(rhoE-0.5*rhoplus*q2);
cplus = sqrt(gamma*p/rhoplus);

vint = uint(2:3)/rhoplus;
rhoEint = uint(4);
vb = vint - dot(vint, n)*n;
qb = norm(vb);
term = 0.5*rhoplus*qb^2;
pb = (gamma-1)*(rhoEint-term);
Fb = [0; pb*n(1); pb*n(2); 0];

udot = dot(vint',n);
s = abs(udot) + cplus;

end