clc
clear
close all

%% =====================================================
% PARAMETERS
%% =====================================================

L  = 1.0;
H  = 0.5;

Nx = 500;
Ny = 500;

U  = 1.0;
Re = 500;

nu = U*L/Re;

%% =====================================================
% GRID
%% =====================================================

x = linspace(0,L,Nx);
y = linspace(0,H,Ny);

[X,Y] = meshgrid(x,y);

%% =====================================================
% BLASIUS-TYPE BOUNDARY LAYER PROFILE
%% =====================================================

u_field = zeros(Ny,Nx);

for j = 1:Nx

    % THEORETICAL BL THICKNESS

    delta = 5*sqrt(nu*x(j)/U);

    if delta < 1e-5
        delta = 1e-5;
    end

    for i = 1:Ny

        eta = y(i)/delta;

        % SMOOTH BOUNDARY LAYER PROFILE

        u_field(i,j) = U*(2*eta - eta^2);

        % LIMIT TO FREESTREAM

        if eta >= 1
            u_field(i,j) = U;
        end

        if u_field(i,j) < 0
            u_field(i,j) = 0;
        end

    end
end

%% =====================================================
% V VELOCITY
%% =====================================================

v_field = zeros(Ny,Nx);

for i = 2:Ny-1

    for j = 2:Nx-1

        v_field(i,j) = -(u_field(i,j+1)-u_field(i,j-1))/(2*(x(2)-x(1)));

    end
end

%% =====================================================
% STREAMFUNCTION
%% =====================================================

psi = zeros(Ny,Nx);

dy = y(2)-y(1);

for j = 1:Nx

    for i = 2:Ny

        psi(i,j) = psi(i-1,j) + u_field(i,j)*dy;

    end
end

%% =====================================================
% VORTICITY
%% =====================================================

omega = zeros(Ny,Nx);

for i = 2:Ny-1

    for j = 2:Nx-1

        omega(i,j) = -(u_field(i+1,j)-u_field(i-1,j))/(2*dy);

    end
end

%% =====================================================
% NUMERICAL BOUNDARY LAYER THICKNESS
%% =====================================================

delta_num = zeros(1,Nx);

for j = 1:Nx

    idx = find(u_field(:,j) >= 0.99*U,1);

    if ~isempty(idx)

        delta_num(j) = y(idx);

    else

        delta_num(j) = NaN;

    end
end

%% =====================================================
% THEORETICAL THICKNESS
%% =====================================================

delta_theory = 5*sqrt(nu*x/U);

%% =====================================================
% VELOCITY CONTOUR
%% =====================================================

figure

contourf(X,Y,u_field,40)

shading interp

colorbar

xlabel('x')
ylabel('y')

title('Velocity Contour')

%% =====================================================
% VORTICITY CONTOUR
%% =====================================================

figure

contourf(X,Y,omega,40)

shading interp

colorbar

xlabel('x')
ylabel('y')

title('Vorticity Contour')

%% =====================================================
% STREAMFUNCTION CONTOUR
%% =====================================================

figure

contour(X,Y,psi,30)

colorbar

xlabel('x')
ylabel('y')

title('Streamfunction Contour')

%% =====================================================
% VELOCITY VECTOR
%% =====================================================

figure

quiver(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),u_field(1:4:end,1:4:end),v_field(1:4:end,1:4:end))

xlabel('x')
ylabel('y')

title('Velocity Vector Plot')

%% =====================================================
% BOUNDARY LAYER THICKNESS
%% =====================================================

figure

plot(x,delta_num,'b','LineWidth',3)

hold on

plot(x,delta_theory,'r--','LineWidth',3)

xlabel('x')

ylabel('\delta')

legend('Numerical','Theory')

title('Boundary Layer Thickness')

grid on

%% =====================================================
% VELOCITY PROFILES
%% =====================================================

figure

hold on
grid on

xloc = [0.1 0.3 0.5 0.7 0.9];

for k = 1:length(xloc)

    [~,j] = min(abs(x-xloc(k)));

    plot(u_field(:,j),y,'LineWidth',2)

end

xlabel('u/U')

ylabel('y')

title('Velocity Profiles')

legend('0.1L','0.3L','0.5L','0.7L','0.9L')