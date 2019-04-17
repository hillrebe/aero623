function [Fb, s] = outflux(uplus, n)

pb = 1;
gamma = 1.4;

rhoplus = uplus(1);
rhou = uplus(2);
rhov = uplus(3);
rhoE = uplus(4);

utemp = rhou/rhoplus;
vtemp = rhov/rhoplus;
q2 = utemp^2 + vtemp^2;
pplus = (gamma-1)*(rhoE-0.5*rhoplus*q2);
Etemp = rhoE/rhoplus;
cplus = sqrt(gamma*pplus/rhoplus);

unplus = dot([utemp vtemp],n);
Jplus = unplus + (2*cplus)/(gamma-1);

Splus = pplus/(rhoplus^gamma);

rhob = (pb/Splus)^(1/gamma);
cb = sqrt(gamma*pb/rhob);

unb = Jplus - 2*cb/(gamma-1);
vbvec = [utemp vtemp] - unplus*n + unb*n;
ub = vbvec(1);
vb = vbvec(2);
qb = norm(vbvec);

rhoEb = pb/(gamma-1) + 0.5*rhob*qb^2;
Hb = rhoEb/rhob + pb/rhob;

Fu_x = [rhob*ub; rhob*ub^2 + pb; rhob*ub*vb; rhob*ub*Hb];
Fu_y = [rhob*vb; rhob*ub*vb; rhob*vb^2 + pb; rhob*vb*Hb];

for i = 1:4
    Fb(i,1) = dot([Fu_x(i) Fu_y(i)], n);
end

s = abs(unplus) + cplus;

end