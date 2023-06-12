function [rho,ux,uy,uz] = push(rho,ux,uy,uz,grid)

% Push with KEP scheme (rho*u*u constraint)
% Forward Euler Timestep
% CD Space derivatives (Unstable?)
% Should be horendously diffusive but conserve KE

%SSP-RK3 (3 stage)
rho0 = rho;
rho_u = zeros(3,grid.Nx,grid.Ny);
rho_u(1,:,:) = rho.*ux;
rho_u(2,:,:) = rho.*uy;
rho_u(3,:,:) = rho.*uz;
rho0_u0 = rho_u;
%Stage 1:
[rho, rho_u] = stage1(rho, rho_u, grid);
%Stage 2:
[rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid);
%Stage 3:
[rho, rho_u] = stage3(rho, rho_u, rho0, rho0_u0, grid);

%Retrieve u:
ux = squeeze(rho_u(1,:,:))./rho;
uy = squeeze(rho_u(2,:,:))./rho;
uz = squeeze(rho_u(3,:,:))./rho;
end


%SSP-RK3 (stage 3)
function [rho, rho_u] = stage3(rho, rho_u, rho0, rho0_u0, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = (1/3)*rho0 + (2/3)*rho_star;
rho_u = (1/3)*rho0_u0 + (2/3)*rho_u_star;

end

%SSP-RK3 (stage 2)
function [rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = (3/4)*rho0 + (1/4)*rho_star;
rho_u = (3/4)*rho0_u0 + (1/4)*rho_u_star;

end

%SSP-RK3 (stage 1)
function [rho, rho_u] = stage1(rho, rho_u, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = rho_star;
rho_u = rho_u_star;

end

function [rho, rho_u] = euler(rho, rho_u, grid)

%Simplify the equations;
c = grid.dt/grid.dx;
R = grid.R;
L = grid.L;
ux = squeeze(rho_u(1,:,:))./rho;
uy = squeeze(rho_u(2,:,:))./rho;
uz = squeeze(rho_u(3,:,:))./rho;
gamma = sqrt(1 + ux.*ux + uy.*uy + uz.*uz);
vx = ux./gamma;
vy = uy./gamma;
vz = uz./gamma;

%Central averages to get ...
puvx_j_plus_half = zeros(3,grid.Nx,grid.Ny);
puvx_j_minus_half = zeros(3,grid.Nx,grid.Ny);
puvy_k_plus_half = zeros(3,grid.Nx,grid.Ny);
puvy_k_minus_half = zeros(3,grid.Nx,grid.Ny);
puvy_j_minus_k_plus_half = zeros(3,grid.Nx,grid.Ny);
puvx_j_plus_half(1,:,:) = 0.5*(rho(R,:).*vx(R,:).*ux(R,:) + rho.*vx.*ux);
puvx_j_plus_half(2,:,:) = 0.5*(rho(R,:).*vx(R,:).*uy(R,:) + rho.*vx.*uy);
puvx_j_plus_half(3,:,:) = 0.5*(rho(R,:).*vx(R,:).*uz(R,:) + rho.*vx.*uz);
puvx_j_minus_half(1,:,:) = 0.5*(rho(L,:).*vx(L,:).*ux(L,:) + rho.*vx.*ux);
puvx_j_minus_half(2,:,:) = 0.5*(rho(L,:).*vx(L,:).*uy(L,:) + rho.*vx.*uy);
puvx_j_minus_half(3,:,:) = 0.5*(rho(L,:).*vx(L,:).*uz(L,:) + rho.*vx.*uz);
puvy_k_plus_half(1,:,:) = 0.5*(rho(:,R).*vy(:,R).*ux(:,R) + rho.*vy.*ux);
puvy_k_plus_half(2,:,:) = 0.5*(rho(:,R).*vy(:,R).*uy(:,R) + rho.*vy.*uy);
puvy_k_plus_half(3,:,:) = 0.5*(rho(:,R).*vy(:,R).*uz(:,R) + rho.*vy.*uz);
puvy_k_minus_half(1,:,:) = 0.5*(rho(:,L).*vy(:,L).*ux(:,L) + rho.*vy.*ux);
puvy_k_minus_half(2,:,:) = 0.5*(rho(:,L).*vy(:,L).*uy(:,L) + rho.*vy.*uy);
puvy_k_minus_half(3,:,:) = 0.5*(rho(:,L).*vy(:,L).*uz(:,L) + rho.*vy.*uz);
puvy_j_minus_k_plus_half(1,:,:) = 0.5*(rho(L,R).*vy(L,R).*ux(L,R) + rho(L,:).*vy(L,:).*ux(L,:));
puvy_j_minus_k_plus_half(2,:,:) = 0.5*(rho(L,R).*vy(L,R).*uy(L,R) + rho(L,:).*vy(L,:).*uy(L,:));
puvy_j_minus_k_plus_half(3,:,:) = 0.5*(rho(L,R).*vy(L,R).*uz(L,R) + rho(L,:).*vy(L,:).*uz(L,:));

%Coef + or -
dv_j_plus = zeros(3,grid.Nx,grid.Ny);
dv_j_minus = zeros(3,grid.Nx,grid.Ny);
dv_k_plus = zeros(3,grid.Nx,grid.Ny);
dv_j_minus_k_plus = zeros(3,grid.Nx,grid.Ny);
dv_j_plus(1,:,:) = (vx - vx(R,:));
dv_j_plus(2,:,:) = (vy - vy(R,:));
dv_j_plus(3,:,:) = (vz - vz(R,:));
dv_k_plus(1,:,:) = (vx - vx(:,R));
dv_k_plus(2,:,:) = (vy - vy(:,R));
dv_k_plus(3,:,:) = (vz - vz(:,R));
dv_j_minus(1,:,:) = (vx(L,:) - vx);
dv_j_minus(2,:,:) = (vy(L,:) - vy);
dv_j_minus(3,:,:) = (vz(L,:) - vz);
dv_j_minus_k_plus(1,:,:) = (vx(L,:) - vx(L,R));
dv_j_minus_k_plus(2,:,:) = (vy(L,:) - vy(L,R));
dv_j_minus_k_plus(3,:,:) = (vz(L,:) - vz(L,R));

%Assemble the coefficients
coef_T2 = zeros(3,grid.Nx,grid.Ny);
coef_T3 = zeros(3,grid.Nx,grid.Ny);
coef_T5 = zeros(3,grid.Nx,grid.Ny);
coef_T6 = zeros(3,grid.Nx,grid.Ny);
coef_T1 = ( (gamma(R,:).*gamma)./(gamma(R,:)-gamma) ).*(1./gamma -  1./gamma(:,R));
coef_T4 = ( (gamma.*gamma(L,:))./(gamma-gamma(L,:)) ).*(1./gamma(L,:) -  1./gamma(L,R));
for i = 1:3
coef_T2(i,:,:) = ( (gamma(R,:).*gamma)./(gamma(R,:)-gamma) ).*squeeze(dv_j_plus(i,:,:));
coef_T3(i,:,:) = ( (gamma(R,:).*gamma)./(gamma(R,:)-gamma) ).*squeeze(dv_k_plus(i,:,:));
coef_T5(i,:,:) = ( (gamma.*gamma(L,:))./(gamma-gamma(L,:)) ).*squeeze(dv_j_minus(i,:,:));
coef_T6(i,:,:) = ( (gamma.*gamma(L,:))./(gamma-gamma(L,:)) ).*squeeze(dv_j_minus_k_plus(i,:,:));
end

%Average pvy
pvy_k_plus_half = 0.5*(rho(:,R).*vy(:,R) + rho.*vy);
pvy_k_minus_half = 0.5*(rho(:,L).*vy(:,L) + rho.*vy);
pvy_j_minus_k__plus_half = 0.5*(rho(L,R).*vy(L,R) + rho(L,:).*vy(L,:));


%Calc terms:
T1 = -coef_T1.*pvy_k_plus_half;
T2 = -squeeze(sum(coef_T2.*puvx_j_plus_half,1));
T3 = -squeeze(sum(coef_T3.*puvy_k_plus_half,1));
T4 = coef_T4.*pvy_j_minus_k__plus_half;
T5 = squeeze(sum(coef_T5.*puvx_j_minus_half,1));
T6 = squeeze(sum(coef_T6.*puvy_j_minus_k_plus_half,1));


%Update rho, u
% Euler Update:
rho = rho - c*(T1 + T2 + T3 + T4 + T5 + T6 + pvy_k_plus_half - pvy_k_minus_half);
rho_u = rho_u - c*(puvx_j_plus_half - puvx_j_minus_half + puvy_k_plus_half - puvy_k_minus_half);

end