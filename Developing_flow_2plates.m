
% DEVELOPING FLOW BETWEEN TWO PARALLEL PLATES



% PARAMETERS


L  = 10;
h  = 1;

U  = 1;

Nx = 50;
Ny = 50;

max_iter = 5000;
tol = 1e-6;


% REYNOLDS NUMBERS

Re_cases = [100 500 1000];


% ENTRY LENGTH STORAGE


Le_numerical = zeros(size(Re_cases));

% LOOP OVER REYNOLDS NUMBERS


for rr = 1:length(Re_cases)

    Re = Re_cases(rr);

    %% GRID

    dx = L/(Nx-1);
    dy = h/(Ny-1);

    x = linspace(0,L,Nx);
    y = linspace(0,h,Ny);

    [X,Y] = meshgrid(x,y);

    %% INITIALIZATION

    psi = zeros(Ny,Nx);
    omega = zeros(Ny,Nx);

    u = zeros(Ny,Nx);
    v = zeros(Ny,Nx);

    %% INLET BC

    u(:,1) = U;

    psi(:,1) = U*y';

    %% RESIDUAL STORAGE

    residual_history = zeros(max_iter,1);

    
    % MAIN ITERATION LOOP
    

    for iter = 1:max_iter

        omega_prev = omega;

        
        % VELOCITY FROM STREAMFUNCTION
       

        u(2:Ny-1,2:Nx-1) = ...
            (psi(3:Ny,2:Nx-1)-psi(1:Ny-2,2:Nx-1))/(2*dy);

        v(2:Ny-1,2:Nx-1) = ...
           -(psi(2:Ny-1,3:Nx)-psi(2:Ny-1,1:Nx-2))/(2*dx);

 
        % CFL TIME STEP
      

        umax = max(abs(u(:))) + 1e-8;
        vmax = max(abs(v(:))) + 1e-8;

        dt_conv = min(dx/umax , dy/vmax);

        dt_diff = Re*min(dx^2,dy^2)/4;

        dt = 0.5*min(dt_conv,dt_diff);

        omega_old = omega;

       
        % INNER ITERATION
        

        for iter2 = 1:300

            psi_old = psi;

            %% UPDATE VELOCITIES

            u(2:Ny-1,2:Nx-1) = ...
                (psi(3:Ny,2:Nx-1)-psi(1:Ny-2,2:Nx-1))/(2*dy);

            v(2:Ny-1,2:Nx-1) = ...
               -(psi(2:Ny-1,3:Nx)-psi(2:Ny-1,1:Nx-2))/(2*dx);

            % VORTICITY TRANSPORT
           

            for i = 2:Ny-1

                for j = 2:Nx-1

                    %% UPWIND X

                    if u(i,j) > 0

                        domega_dx = ...
                            (omega_old(i,j)-omega_old(i,j-1))/dx;

                    else

                        domega_dx = ...
                            (omega_old(i,j+1)-omega_old(i,j))/dx;

                    end

                    %% UPWIND Y

                    if v(i,j) > 0

                        domega_dy = ...
                            (omega_old(i,j)-omega_old(i-1,j))/dy;

                    else

                        domega_dy = ...
                            (omega_old(i+1,j)-omega_old(i,j))/dy;

                    end

                    %% CONVECTION

                    conv = ...
                        u(i,j)*domega_dx + ...
                        v(i,j)*domega_dy;

                    %% DIFFUSION

                    lap_omega = ...
                        (omega_old(i+1,j)-2*omega_old(i,j)+omega_old(i-1,j))/dy^2 ...
                        + ...
                        (omega_old(i,j+1)-2*omega_old(i,j)+omega_old(i,j-1))/dx^2;

                    diff = lap_omega/Re;

                    %% TIME UPDATE

                    omega(i,j) = ...
                        omega_old(i,j) + dt*(diff-conv);

                end
            end

         
            % BOUNDARY CONDITIONS
            
            % TOP WALL

            omega(end,:) = ...
                2*(U*h - psi(end-1,:))/dy^2;

            % SYMMETRY LINE

            omega(1,:) = 0;

            % INLET

            omega(:,1) = 0;

            % OUTLET

            omega(:,end) = omega_old(:,end-1);

            %% PSI BC

            psi(1,:) = 0;

            psi(end,:) = U*h;

            psi(:,1) = U*y';

            psi(:,end) = psi(:,end-1);

            %% =========================================
            % POISSON EQUATION
            %% =========================================

            for i = 2:Ny-1

                for j = 2:Nx-1

                    psi(i,j) = ((psi(i+1,j)+psi(i-1,j))/dy^2 + (psi(i,j+1)+psi(i,j-1))/dx^2 + omega(i,j)) / (2/dx^2 + 2/dy^2);

                end
            end

            %% INNER CONVERGENCE

            if max(abs(psi(:)-psi_old(:))) < tol
                break;
            end

        end

        %% =============================================
        % GLOBAL CONVERGENCE
        %% =============================================

        res = max(abs(omega(:)-omega_prev(:)));

        residual_history(iter) = res;

        if mod(iter,100)==0

            fprintf( ...
            'Re = %d  Iter = %d  Residual = %.3e\n', ...
            Re,iter,res);

        end

        if res < tol

            fprintf( ...
            '\nSteady state reached for Re = %d at iteration %d\n', ...
            Re,iter);

            break;

        end

    end

   
    % FINAL VELOCITY


    u(2:Ny-1,2:Nx-1) = (psi(3:Ny,2:Nx-1)-psi(1:Ny-2,2:Nx-1))/(2*dy);


    % ENTRY LENGTH
   

    u_centerline = u(round(Ny/2),:);

    u_fully = max(u(:,end));

    idx = find(u_centerline >= 0.99*u_fully,1);

    if ~isempty(idx)

        Le_numerical(rr) = x(idx);

    else

        Le_numerical(rr) = NaN;

    end

    
    % VELOCITY PROFILES
   

    figure;

    hold on;
    grid on;

    x_locs = [0.05 0.2 0.95]*L;

    colors = lines(length(x_locs));

    for k = 1:length(x_locs)

        [~,j] = min(abs(x-x_locs(k)));

        plot( ...
            u(:,j), ...
            y, ...
            'LineWidth',2, ...
            'Color',colors(k,:), ...
            'DisplayName', ...
            sprintf('x = %.2f',x(j)));

    end

    %% ANALYTICAL PROFILE

    y_anal = linspace(0,h,200);

    u_anal = 1.5*U*(1-(y_anal/h).^2);

    plot( ...
        u_anal, ...
        y_anal, ...
        'k--', ...
        'LineWidth',2, ...
        'DisplayName','Analytical');

    xlabel('Velocity');

    ylabel('y');

    title(sprintf('Velocity Profiles (Re = %d)',Re));

    legend show;

   
    % BOUNDARY LAYER THICKNESS
   

    delta_num = zeros(1,Nx);

    for j = 1:Nx

        target = 0.99*max(u(:,end));

        idx2 = find(u(:,j) >= target,1);

        if ~isempty(idx2)

            delta_num(j) = y(idx2);

        else

            delta_num(j) = NaN;

        end

    end

    %% THEORETICAL THICKNESS

    delta_theory = 5*sqrt(x/Re);

   
    % BOUNDARY LAYER THICKNESS PLOT

    figure

    plot(x,delta_num,'b','LineWidth',3)

    hold on

    plot(x,delta_theory,'r--','LineWidth',3)

    xlabel('x')

    ylabel('\delta')

    legend('Numerical','Theory')

    title(sprintf( ...
        'Boundary Layer Thickness (Re = %d)', ...
        Re))

    grid on

   
    % RESIDUAL PLOT
    

    figure;

    semilogy( ...
        1:iter, ...
        residual_history(1:iter), ...
        'LineWidth',2);

    xlabel('Iteration');

    ylabel('Residual');

    title(sprintf('Residual vs Iteration (Re = %d)',Re));

    grid on;

end


% ENTRY LENGTH VS RE


Le_theory = 0.05*Re_cases*(2*h);

figure;

plot( ...
    Re_cases, ...
    Le_numerical, ...
    'bo-','LineWidth',2);

hold on;

plot( ...
    Re_cases, ...
    Le_theory, ...
    'r--','LineWidth',2);

xlabel('Reynolds Number');

ylabel('Entry Length');

title('Entry Length vs Reynolds Number');

legend('Numerical','Theory');

grid on;