function [u_final, Rnorm] = FVM(meshdata, Nmax, order, boundaries)

% boundaries = 0: freestream
% boundaries = 1: use boundaries

E2N = meshdata.E2N;
% V = meshdata.V;
% Area = meshdata.area;
% I2E = meshdata.I2E;
% B2E = meshdata.B2E;
% In = meshdata.In;
% Bn = meshdata.Bn;

% Setup
tol = 1e-7;
R = 1;
p_inf = 1;
M_inf = 0.2;
gamma = 1.4;
T_t = 1 + (gamma-1)/2*M_inf^2;
% p_t = T_t^(gamma/(gamma-1));
cp = (gamma/(gamma-1))*R;
T = T_t/(1 + (gamma-1)/2*M_inf^2);
h = cp*T;
c = sqrt((gamma-1)*h);
v_IC = 0;
u_IC = M_inf*c;
rho_inf = gamma*p_inf/c^2;
q_IC = sqrt(u_IC^2+v_IC^2);
E = (p_inf/(gamma-1)+0.5*rho_inf*q_IC^2)/rho_inf;
% H = E+p_inf/rho_inf;

u_vec = [rho_inf; rho_inf*u_IC; rho_inf*v_IC; rho_inf*E];

nelem = length(E2N);
% nInt = length(I2E);
% nBound = length(B2E);

Rnorm = zeros(Nmax,1);

% RK2 Freestream
if boundaries == 0
    
    umat = ones(nelem,4);
    umat(:,1) = umat(:,1)*u_vec(1);
    umat(:,2) = umat(:,2)*u_vec(2);
    umat(:,3) = umat(:,3)*u_vec(3);
    umat(:,4) = umat(:,4)*u_vec(4);
    
    if order == 1
        for n = 1:Nmax
            [Run, dtA] = residual_1(meshdata, umat, u_vec, boundaries);
            uiFE = umat - dtA.*Run;
%             [RuFE, ~] = residual_1(meshdata, uiFE, u_vec, boundaries);
%             un1 = 0.5*(umat + uiFE - dtA.*RuFE);
%             umat = un1;
            umat = uiFE;
            Rnorm(n) = max(max(abs(Run)));
        end
    elseif order == 2
        for n = 1:Nmax
            [Run, dtA] = residual_2(meshdata, umat, u_vec, boundaries);
            uiFE = umat - dtA.*Run;
%             [RuFE, ~] = residual_2(meshdata, uiFE, u_vec, boundaries);
%             un1 = 0.5*(umat + uiFE - dtA.*RuFE);
%             umat = un1;
            umat = uiFE;
            Rnorm(n) = max(max(abs((Run))));
        end  
    end

% RK2 with boundaries
elseif boundaries == 1
    
%     umat = ones(nelem,4);
%     umat(:,1) = umat(:,1)*u_vec(1);
%     umat(:,2) = umat(:,2)*u_vec(2);
%     umat(:,3) = umat(:,3)*u_vec(3);
%     umat(:,4) = umat(:,4)*u_vec(4);
    load('u_2412.mat');
    umat = u_2412;

    if order == 1
        for n = 1:Nmax
            [Run, dtA] = residual_1(meshdata, umat, u_vec, boundaries);
            uiFE = umat - dtA.*Run;
%             [RuFE, ~] = residual_1(meshdata, uiFE, u_vec, boundaries);
%             un1 = 0.5*(umat + uiFE - dtA.*RuFE);
%             umat = un1;
            umat = uiFE;
            Rnorm(n) = max(max(abs(Run)));
            
            if floor(n/100) == n/100
                save('utemp','umat');
                fprintf('Iterations: %g \n',n);
                fprintf('Residual: %g \n \n',Rnorm(n));
            end
            
            if Rnorm(n) < tol
                Rnorm(n+1:end) = [];
                save('ufinal','umat');
                break;
            end
        end
    elseif order == 2
        for n = 1:Nmax
            [Run, dtA] = residual_2(meshdata, umat, u_vec, boundaries);
            uiFE = umat - dtA.*Run;
%             [RuFE, ~] = residual_2(meshdata, uiFE, u_vec, boundaries);
%             un1 = 0.5*(umat + uiFE - dtA.*RuFE);
%             umat = un1;
            umat = uiFE;
            Rnorm(n) = max(max(abs(Run)));
            
            if floor(n/100) == n/100
                save('utemp','umat');
                fprintf('Iterations: %g \n',n);
                fprintf('Residual: %g \n \n',Rnorm(n));
            end
            
            if Rnorm(n) < tol
                Rnorm(n+1:end) = [];
                save('ufinal','umat');
                break;
            end
        end
    end
    
else
    fprintf('Invalid boundaries input');
end

u_final = umat;
% Rhist = Rnorm;

end