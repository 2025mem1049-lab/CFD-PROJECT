clear; clc; close all;

%% Parameters
L = 1;
Re = 100;
Nx = 51;
Ny = 51;
U = 1;

max_iter = 5000;
tol = 1e-6;

%% Grid
dx = L/(Nx-1);
dy = L/(Ny-1);

x = linspace(0,L,Nx);
y = linspace(0,L,Ny);

[X,Y] = meshgrid(x,y);

%% Initialization
psi = zeros(Ny,Nx);
omega = zeros(Ny,Nx);

u = zeros(Ny,Nx);
v = zeros(Ny,Nx);

%% Residual Storage
residual_history = zeros(max_iter,1);

%% Time Loop
for iter = 1:max_iter

    omega_old = omega;
    psi_old_global = psi;

    %% Velocity from streamfunction
    u(2:end-1,2:end-1) = ...
        (psi(3:end,2:end-1) - psi(1:end-2,2:end-1))/(2*dy);

    v(2:end-1,2:end-1) = ...
       -(psi(2:end-1,3:end) - psi(2:end-1,1:end-2))/(2*dx);

    %% Stable time step
    umax = max(abs(u(:))) + 1e-8;
    vmax = max(abs(v(:))) + 1e-8;

    dt_conv = min(dx/umax , dy/vmax);

    dt_diff = Re*dx^2/4;

    dt = 0.5*min(dt_conv,dt_diff);

    %% --- Vorticity Transport Equation ---
    for i = 2:Ny-1

        for j = 2:Nx-1

            % Convection terms
            domega_dx = ...
                (omega_old(i,j+1)-omega_old(i,j-1))/(2*dx);

            domega_dy = ...
                (omega_old(i+1,j)-omega_old(i-1,j))/(2*dy);

            conv = ...
                u(i,j)*domega_dx + ...
                v(i,j)*domega_dy;

            % Diffusion term
            lap_omega = ...
                (omega_old(i+1,j)-2*omega_old(i,j)+omega_old(i-1,j))/dy^2 + ...
                (omega_old(i,j+1)-2*omega_old(i,j)+omega_old(i,j-1))/dx^2;

            diff = lap_omega/Re;

            % Time update
            omega(i,j) = ...
                omega_old(i,j) + dt*(diff-conv);

        end
    end

    %% ===================================
    %% Boundary Conditions
    %% ===================================

    % Streamfunction BC
    psi(1,:) = 0;
    psi(end,:) = 0;
    psi(:,1) = 0;
    psi(:,end) = 0;

    %% Vorticity BC

    % Bottom wall
    omega(1,:) = -2*psi(2,:)/dy^2;

    % Top moving lid
    omega(end,:) = ...
        -2*psi(end-1,:)/dy^2 ...
        -2*U/dy;

    % Left wall
    omega(:,1) = -2*psi(:,2)/dx^2;

    % Right wall
    omega(:,end) = -2*psi(:,end-1)/dx^2;

    %% --- Solve Poisson Equation ---
    for k = 1:200

        psi_old = psi;

        for i = 2:Ny-1

            for j = 2:Nx-1

                psi(i,j) = ...
                    ((psi(i+1,j)+psi(i-1,j))/dy^2 + ...
                    (psi(i,j+1)+psi(i,j-1))/dx^2 + ...
                    omega(i,j)) ...
                    /(2/dx^2 + 2/dy^2);

            end
        end

        if max(abs(psi(:)-psi_old(:))) < tol
            break;
        end
    end

    %% Convergence Check
    res = max(abs(omega(:)-omega_old(:)));

    residual_history(iter) = res;

    if mod(iter,100)==0
        fprintf('Iter = %d   Residual = %.3e\n',iter,res);
    end

    if res < tol
        fprintf('\nSteady state reached\n');
        break;
    end

end

%% Final velocity
u(2:end-1,2:end-1) = ...
    (psi(3:end,2:end-1)-psi(1:end-2,2:end-1))/(2*dy);

v(2:end-1,2:end-1) = ...
   -(psi(2:end-1,3:end)-psi(2:end-1,1:end-2))/(2*dx);

%% =====================================================
% STREAM FUNCTION CONTOUR
%% =====================================================

figure;

contourf(X,Y,psi,30,'LineStyle','none');

colorbar;

title('Stream Function');

xlabel('x');
ylabel('y');

%% =====================================================
% VORTICITY CONTOUR
%% =====================================================

figure;

contourf(X,Y,omega,30,'LineStyle','none');

colorbar;

title('Vorticity Contour');

xlabel('x');
ylabel('y');

%% =====================================================
% VELOCITY VECTOR
%% =====================================================

figure;

quiver( ...
    X(1:2:end,1:2:end), ...
    Y(1:2:end,1:2:end), ...
    u(1:2:end,1:2:end), ...
    v(1:2:end,1:2:end));

title('Velocity Field');

xlabel('x');
ylabel('y');

%% =====================================================
% U VELOCITY CONTOUR
%% =====================================================

figure;

contourf(X,Y,u,30,'LineStyle','none');

colorbar;

title('U Velocity Contour');

xlabel('x');
ylabel('y');

%% =====================================================
% RESIDUAL VS ITERATION
%% =====================================================

figure;

semilogy( ...
    1:iter, ...
    residual_history(1:iter), ...
    'LineWidth',2);

xlabel('Iteration');

ylabel('Residual');

title('Residual vs Iteration');

grid on;