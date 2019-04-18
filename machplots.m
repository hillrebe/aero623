function machplots(umat, meshdata)

V = meshdata.V;
E2N = meshdata.E2N;

figure

gamma = 1.4;

for i = 1:length(umat)
    rho = umat(i,1);
    rhou = umat(i,2);
    rhov = umat(i,3);
    rhoE = umat(i,4);

    utemp = rhou/rho;
    vtemp = rhov/rho;
    q2 = utemp^2 + vtemp^2;
    p(i) = (gamma-1)*(rhoE-0.5*rho*q2);
    c = sqrt(gamma*p(i)/rho);
    M(i) = sqrt(q2)/c;
    
end

Mmin = min(M);
Mmax = max(M);
colorref = parula(1000);

for i = 1:length(umat)
    nodes = E2N(i,:);
%     xtot = 0;
%     ytot = 0;
    
    for j = 1:4
%         xtot = xtot + V(nodes(j),1);
%         ytot = ytot + V(nodes(j),2);
        xvec(j) = V(nodes(j),1);
        yvec(j) = V(nodes(j),2);
    end
    
%     cent_x(i) = xtot/4;
%     cent_y(i) = ytot/4;
    
%     frac = 1000*(M(i)-Mmin)/(Mmax-Mmin);
%     frac_int = round(frac);
%     if frac_int == 0
%         frac_int = 1;
%     end
%     color = colorref(frac_int,:);
    
    fill(xvec,yvec,M(i),'edgecolor','none')
    hold on
end

axis image
bar = colorbar;
caxis([Mmin, Mmax])
xlabel('x-position')
ylabel('y-position')
h = colorbar;
ylabel(h,'Mach Number')

figure

for i = 1:length(umat)
    nodes = E2N(i,:);

    for j = 1:4
        xvec(j) = V(nodes(j),1);
        yvec(j) = V(nodes(j),2);
    end
    
    fill(xvec,yvec,p(i),'edgecolor','none')
    hold on
end

axis image
h = colorbar;
ylabel(h,'Pressure')
% bar = colorbar;
% caxis([Mmin, Mmax])
xlabel('x-position')
ylabel('y-position')

end