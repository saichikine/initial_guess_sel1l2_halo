function results = multiple_shooting(free_vars_guess, constraint_handle, jacobian_handle, varargin)
    
    %% Input handling
    
    default_tol = 1e-12;
    default_max_iter = 30;
%     default_bool_plot = 0;
    default_bool_verbose = 1;
    default_bool_timer = 1;
    
    p = inputParser;
    valid_scalar_pos_num = @(x) isnumeric(x) && isscalar(x) && (x>0);
    valid_ICs = @(x) all(size(initial_guess_ICs==[6 1]));
    valid_f_handle = @(x) isa(x, 'function_handle');
    valid_bool = @(x) isa(x, 'logical');
    valid_int = @(x) mod(x,1)==0;
    
    addRequired(p,'free_vars_guess',@isnumeric);
    addRequired(p,'constraint_handle',valid_f_handle);
    
%     addParameter(p,'bool_plot',default_bool_plot,valid_bool);
%     addParameter(p,'plot_handle', valid_f_handle);
    addParameter(p,'bool_verbose',default_bool_verbose,valid_bool);
    addParameter(p,'bool_timer',default_bool_timer,valid_bool);
    addParameter(p,'tol',default_tol,valid_scalar_pos_num);
    addParameter(p,'max_iter',default_max_iter,valid_int);
    
    parse(p,free_vars_guess,constraint_handle,jacobian_handle,varargin{:});
    
    free_vars_guess = p.Results.free_vars_guess;
    constraint_handle = p.Results.contstraints_handle;
    jacobian_handle = p.Results.jacobian_handle;
    
%     bool_plot = p.Results.bool_plot;
    bool_verbose = p.Results.bool_verbose;
    bool_timer = p.Results.bool_timer;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
 
    %% Setup

    % Build initial big X
    chi = free_vars_guess;
    
    % Initial constraint vector, chi is free variable vector
    g = constraint_handle(chi);

    %% Multiple Shooting Loop
    
    % Stop condition for loop
    bool_converged = false;
    
    gamma_vec = linspace(0,1.5,20);
    
    err_mag = 100;
    err_mag_old = err_mag;
    err_mag_hist = [err_mag];

    counter = 0;
    
    if bool_verbose
        fprintf("Starting Multiple Shooter.\n")
    end
    if bool_timer
        tic
    end
    while ~bool_converged
        
        % If verbose messages are desired, print info
        if bool_verbose
            fprintf("Iteration %i\nCurrent residual magnitude is: %d\n",counter,err_mag);
        end
        
        % Compute Jacobian
        dgdchi = jacobian_handle(chi);
        
        if det(dgdchi==0)
            error("Jacobian not invertible.")
        end

        % Test updates with different values of gamma to find greatest error reduction
        err_mag_vec = NaN(1,length(gamma_vec));
        gtest = g;
        parfor c = 1:length(gamma_vec)
            delta_chi = dgdchi'/(dgdchi*dgdchi')*(-gamma_vec(c)*g);
            chi_test = free_vars + delta_chi;

            gtest = constraint_handle(chi_test);
            err_mag = norm(gtest);
            err_mag_vec(c) = err_mag;
        end
        
        % Find correction that reduces error the most
        min_gamma_index = find(err_mag_vec == min(err_mag_vec));
        err_mag = err_mag_vec(min_gamma_index);
        
        % If error did not decrease, stop loop
        if err_mag >= err_mag_old
            err_mag_hist = [err_mag_hist, err_mag_old];
            fprintf("Residual did not decrease. Terminating with a final residual of %d.\n", err_mag_old)
            break
        else
            err_mag_old = err_mag; % If error continues to decrease, reset old error for comparison
        end
        
        % Use gamma corresponding to maximum error decrease
        gamma = gamma_vec(min_gamma_index);

        % Compute update to free variable vector chi
        delta_chi = dgdchi'/(dgdchi*dgdchi')*(-gamma*g);
        
        % Update free variables
        chi = chi + delta_chi;
        
        % Compute constraint vector
        g = constraint_handle(chi);
        
        % Compute error
        err_mag = norm(g);
        
        % Save error magnitude
        err_mag_hist = [err_mag_hist, err_mag];

        if err_mag <= tol
            bool_converged = 1;
            fprintf("Converged, with a final error of %d\n", err_mag);
        elseif counter > max_iter
            bool_converged = 1;
            fprintf("Exceeded maximum number of iterations (%d)\n", max_iter);
        end
    end
    
    % If timer, stop timer
    if bool_timer
        toc
    end

    %% Plot comparison (starting off)
    
    if bool_plot
        
        L_points = lagrangePoints(mu);
        x_L1 = L_points(1,1);
        x_L2 = L_points(1,2);

        %% Plot comparison (after done)
        
        % 3D Plot
%         figure
%         addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
%         plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
%         plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
%         plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
%         plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
%         arcs = {};
%         for i = 1:N-1
%             arcs{i} = scatter3(X_hists{i}(:,1), X_hists{i}(:,2), X_hists{i}(:,3),'.'); hold on
%         end
%         title('Multiple Shooting Converged')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         grid on;
%         legend();
%         hold off
        
        % 2D Plot
        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot(1-mu, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot(X_hist_total(1,1), X_hist_total(1,2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot(x_L1, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot(x_L2, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        arcs = {};
        for i = 1:N-1
            arcs{i} = scatter(X_hists{i}(:,1), X_hists{i}(:,2), '.'); hold on
        end
        title('Multiple Shooting Converged')
        xlabel('X')
        ylabel('Y')
        axis('equal')
        grid on;
        legend();
        hold off

%         figure
%         addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
%         plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
%         plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
%         plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
%         plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
%         plot3(X_hist_total(:,1), X_hist_total(:,2), X_hist_total(:,3), 'r-','DisplayName', 'Converged Connection'); hold on
%         title('Multiple Shooting Converged')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         grid on;
%         legend();
%         hold off
    end
    
    %% Save results
    
    results = struct('X0',orbit_IC,'total_time',total_time,'jacobi_constant',C,'free_vars',chi);
end