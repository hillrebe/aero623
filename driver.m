clear all
close all

global E2N I2E B2E In Area V Bn

%% Setup
CFL_FE = 0.9;
tol = 1e-7;
Nmax = 10000;
R = 1;
p_inf = 1;
M_inf = 0.5;
gamma = 1.4;
T_t = 1 + (gamma-1)/2*M_inf^2;
p_t = T_t^(gamma/(gamma-1));
alpha = 0;
cp = (gamma/(gamma-1))*R;
T = T_t/(1 + (gamma-1)/2*M_inf^2);
h = cp*T;
c = sqrt((gamma-1)*h);
v_IC = 0;
u_IC = M_inf*c;
rho_inf = gamma*p_inf/c^2;
q_IC = sqrt(u_IC^2+v_IC^2);
E = (p_inf/(gamma-1)+0.5*rho_inf*q_IC^2)/rho_inf;
H = E+p_inf/rho_inf;

kmax = 1;

cl = zeros(kmax,2);
cd = zeros(kmax,2);
dof = zeros(kmax,1);
Es = zeros(kmax,2);

for k = 1:kmax
    %% Initialization
    fileOpen = sprintf('bump%g.gri',k-1);
    plotgri(fileOpen);

    u_vec = [rho_inf; rho_inf*u_IC; rho_inf*v_IC; rho_inf*E];

    nelem = length(E2N);
    nInt = length(I2E);
    nBound = length(B2E);
    

    %% FE Freestream

%     umat = ones(nelem,4);
%     umat(:,1) = umat(:,1)*u_vec(1);
%     umat(:,2) = umat(:,2)*u_vec(2);
%     umat(:,3) = umat(:,3)*u_vec(3);
%     umat(:,4) = umat(:,4)*u_vec(4);
%     
%     for n = 1:Nmax
%         
%         Rvec = zeros(nelem,4);
%         si = zeros(nelem,1);
%         
%         for i = 1:nInt
%             eL = I2E(i,1);
%             eR = I2E(i,3);
%             normal = In(i,:);
%             uL = umat(eL,:)';
%             uR = umat(eR,:)';
%             [Fhat, s] = roeflux(uL,uR,normal);
%             
%             faceL = I2E(i,2);
%             faceR = I2E(i,4);
%             nodes = E2N(eL,:);
%             nodes(faceL) = [];
%             dx = V(nodes(2),1) - V(nodes(1),1);
%             dy = V(nodes(2),2) - V(nodes(1),2);
%             deltal(eL,faceL) = sqrt(dx^2+dy^2);
%             deltal(eR,faceR) = sqrt(dx^2+dy^2);
%             smat(eL,faceL) = abs(s);
%             smat(eR,faceR) = abs(s);
%             
%             Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal(eL,faceL);
%             Rvec(eR,:) = Rvec(eR,:) - Fhat'*deltal(eL,faceL);
%         end
%         
%         for i = 1:nBound
%             eL = B2E(i,1);
%             normal = Bn(i,:);
%             uL = umat(eL,:)';
%             uR = u_vec;
%             [Fhat s] = roeflux(uL,uR,normal);
%             
%             faceL = B2E(i,2);
%             nodes = E2N(eL,:);
%             nodes(faceL) = [];
%             dx = V(nodes(2),1) - V(nodes(1),1);
%             dy = V(nodes(2),2) - V(nodes(1),2);
%             deltal(eL,faceL) = sqrt(dx^2+dy^2);
%             smat(eL,faceL) = abs(s);
%             
%             Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal(eL,faceL);
%         end
%         
%         for i = 1:nelem
%             sie = smat(i,:);
%             dlie = deltal(i,:);
%             product = sie.*dlie;
%             si(i,1) = sum(product);
%         end
%         
%         dtA = 2*CFL_FE./si;
%     %     dtA = 0.01;
%         umat = umat - dtA.*Rvec;
%         Rnorm(n) = norm(Rvec);
%         
%     end
%     
%     figure
%     plot(1:Nmax, Rnorm)
%     xlabel('Number of iterations')
%     ylabel('Maximum residual value')
%     Rnorm(1)

    %% FE with boundaries

    umat_FE = ones(nelem,4);
    umat_FE(:,1) = umat_FE(:,1)*u_vec(1);
    umat_FE(:,2) = umat_FE(:,2)*u_vec(2);
    umat_FE(:,3) = umat_FE(:,3)*u_vec(3);
    umat_FE(:,4) = umat_FE(:,4)*u_vec(4);
    
    Rnorm = zeros(Nmax,1);

    for n = 1:Nmax

        Rvec = zeros(nelem,4);
        sdl = zeros(nelem,1);

        for i = 1:nInt
            eL = I2E(i,1);
            eR = I2E(i,3);
            normal = In(i,:);
            uL = umat_FE(eL,:)';
            uR = umat_FE(eR,:)';
            [Fhat, s] = roeflux(uL,uR,normal);

            faceL = I2E(i,2);
            faceR = I2E(i,4);
            nodes = E2N(eL,:);
            nodes(faceL) = [];
            dx = V(nodes(2),1) - V(nodes(1),1);
            dy = V(nodes(2),2) - V(nodes(1),2);
%             deltal(eL,faceL) = sqrt(dx^2+dy^2);
%             deltal(eR,faceR) = sqrt(dx^2+dy^2);
%             smat(eL,faceL) = abs(s);
%             smat(eR,faceR) = abs(s);
            deltal = sqrt(dx^2+dy^2);
            sdl(eL) = sdl(eL) + s*deltal;
            sdl(eR) = sdl(eR) + s*deltal;

            Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal;
            Rvec(eR,:) = Rvec(eR,:) - Fhat'*deltal;
        end

        for i = 1:nBound
            eL = B2E(i,1);
            normal = Bn(i,:);
            uL = umat_FE(eL,:)';
            uR = u_vec;
    %         [Fhat s] = roeflux(uL,uR,normal);

            if abs(normal(2)) > 0
                [Fhat, s] = wallflux(uL, normal);
            elseif normal(1) < 0
                [Fhat, s] = influx(uL, normal);
            else
                [Fhat, s] = outflux(uL, normal);
            end           

            faceL = B2E(i,2);
            nodes = E2N(eL,:);
            nodes(faceL) = [];
            dx = V(nodes(2),1) - V(nodes(1),1);
            dy = V(nodes(2),2) - V(nodes(1),2);
%             deltal(eL,faceL) = sqrt(dx^2+dy^2);
%             smat(eL,faceL) = abs(s);
            deltal = sqrt(dx^2+dy^2);
            sdl(eL) = sdl(eL) + s*deltal;

            Rvec(eL,:) = Rvec(eL,:) + Fhat'*deltal;
        end

%         for j = 1:nelem
%             sie = smat(j,:);
%             dlie = deltal(j,:);
%             product = sie.*dlie;
%             si(j,1) = sum(product);
%         end

        dtA = 2*CFL_FE./sdl;
    %     dtA = 0.01;
        Rnorm(n) = max(max(Rvec));
        umat_FE = umat_FE - dtA.*Rvec;

        if Rnorm(n) < tol
            Rnorm(n+1:end) = [];
            break
        end  

    end

    [cl(k,1), cd(k,1), Es(k,1), ~] = calcOutputs(umat_FE);
    dof(k) = nelem;
    figure
    loglog(1:n, Rnorm)
    xlabel('Number of iterations')
    ylabel('Maximum residual value')
    
    %% RK2 Freestream
    
%     umat = ones(nelem,4);
%     umat(:,1) = umat(:,1)*u_vec(1);
%     umat(:,2) = umat(:,2)*u_vec(2);
%     umat(:,3) = umat(:,3)*u_vec(3);
%     umat(:,4) = umat(:,4)*u_vec(4);
%     
%     for n = 1:Nmax
%         [Run, dtA] = residual(umat, u_vec, 0);
%         uiFE = umat - dtA.*Run;
%         [RuFE, ~] = residual(uiFE, u_vec, 0);
%         un1 = 0.5*(umat + uiFE - dtA.*RuFE);
%         umat = un1;
%         Rnorm(n) = max(max((Run)));
%     end
%     
%     figure
%     plot(1:Nmax, Rnorm)
%     xlabel('Number of iterations')
%     ylabel('Maximum residual value')
%     Rnorm(1)
%     
    
   %% RK2 with boundaries
    
%     umat_RK2 = ones(nelem,4);
%     umat_RK2(:,1) = umat_RK2(:,1)*u_vec(1);
%     umat_RK2(:,2) = umat_RK2(:,2)*u_vec(2);
%     umat_RK2(:,3) = umat_RK2(:,3)*u_vec(3);
%     umat_RK2(:,4) = umat_RK2(:,4)*u_vec(4);
    umat_RK2 = umat_FE;
    
    for n = 1:Nmax
        [Run, dtA] = residual(umat_RK2, u_vec, 1);
        uiFE = umat_RK2 - dtA.*Run;
        [RuFE, ~] = residual(uiFE, u_vec, 1);
        un1 = 0.5*(umat_RK2 + uiFE - dtA.*RuFE);
%         un1 = umat_RK2 - 0.5*dtA.*(Run+RuFE);
        umat_RK2 = un1;
        Rnorm_RK2(n) = max(max(Run));
        if Rnorm_RK2(n) < tol
            Rnorm_RK2(n+1:end) = [];
            break;
        end
    end
    
    [cl(k,2), cd(k,2), Es(k,2), ~] = calcOutputs(umat_RK2);
    
    figure
    loglog(1:length(Rnorm_RK2), Rnorm_RK2)
    xlabel('Number of iterations')
    ylabel('Maximum residual value')
end
% 
%% Calculate errors

% cl_exact = 1.537095;
% cd_exact = 2.94278e-6;
% cl_error = cl - cl_exact;
% cd_error = cd - cd_exact;
% 
% if kmax > 1
%     slope_cl = zeros(kmax-1,2);
%     slope_cd = zeros(kmax-1,2);
%     slope_es = zeros(kmax-1,2);
%     for m = 1:(kmax-1)
%         slope_cl(m,1) = calcSlope(cl_error(m,1),cl_error(m+1,1),dof(m),dof(m+1));
%         slope_cd(m,1) = calcSlope(cd_error(m,1),cd_error(m+1,1),dof(m),dof(m+1));
%         slope_es(m,1) = calcSlope(Es(m,1),Es(m+1,1),dof(m),dof(m+1));
%         slope_cl(m,2) = calcSlope(cl_error(m,2),cl_error(m+1,2),dof(m),dof(m+1));
%         slope_cd(m,2) = calcSlope(cd_error(m,2),cd_error(m+1,2),dof(m),dof(m+1));
%         slope_es(m,2) = calcSlope(Es(m,2),Es(m+1,2),dof(m),dof(m+1));
%     end
% end
% 
% %% Plots
% 
% figure
% loglog(sqrt(dof), abs(cl_error(:,1)),'o-')
% hold on
% loglog(sqrt(dof), abs(cl_error(:,2)),'o-')
% xlabel('$$\sqrt{dof}$$','Interpreter','latex')
% ylabel('Error in c_L')
% legend('Forward Euler','RK2')
% 
% figure
% loglog(sqrt(dof), cd_error(:,1),'o-')
% hold on
% loglog(sqrt(dof), cd_error(:,2),'o-')
% xlabel('$$\sqrt{dof}$$','Interpreter','latex')
% ylabel('Error in c_D')
% legend('Forward Euler','RK2')
% 
% figure
% loglog(sqrt(dof), Es(:,1),'o-')
% hold on
% loglog(sqrt(dof), Es(:,2),'o-')
% xlabel('$$\sqrt{dof}$$','Interpreter','latex')
% ylabel('Error in E_s')
% legend('Forward Euler','RK2')