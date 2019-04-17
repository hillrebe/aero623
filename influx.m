function [Fb, s] = influx(uplus, n)

Minf = 0.5;
gamma = 1.4;
R = 1;
Tt = 1+(gamma-1)/2*Minf^2;
pt = Tt^(gamma/(gamma-1));
alpha = 0;

rho = uplus(1);
rhou = uplus(2);
rhov = uplus(3);
rhoE = uplus(4);

utemp = rhou/rho;
vtemp = rhov/rho;
q2 = utemp^2 + vtemp^2;
p = (gamma-1)*(rhoE-0.5*rho*q2);
Etemp = rhoE/rho;
cplus = sqrt(gamma*p/rho);

unplus = dot([utemp vtemp],n);
Jplus = unplus + (2*cplus)/(gamma-1);

nin = [cos(alpha) sin(alpha)];
dn = dot(nin, n);

term1 = gamma*R*Tt*dn^2 - ((gamma-1)/2)*Jplus^2;
term2 = 4*gamma*R*Tt*dn/(gamma-1);
term3 = 4*gamma*R*Tt/(gamma-1)^2;
term4 = Jplus^2;

Mbvec = roots([term1, term2, term3-term4]);
Mb_poss = Mbvec(Mbvec>=0);
Mb = min(Mb_poss);

Tb = Tt/(1+0.5*(gamma-1)*Mb^2);

pb = pt*(Tb/Tt)^(gamma/(gamma-1));
rhob = pb/(R*Tb);
cb = sqrt(gamma*pb/rhob);
vbvec = Mb*cb*nin;
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