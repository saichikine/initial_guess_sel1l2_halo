function results = shooter(free_vars_guess, constraint_func, varargin)
    
    %% Input handling
    
    default_tol = 1e-12;
    default_max_iter = 30;
    default_bool_verbose = true;
    default_bool_timer = true;
    default_bool_full_hist = false;
    default_jacobian_func = @(chi) false;
    
    p = inputParser;
    valid_scalar_pos_num = @(x) isnumeric(x) && isscalar(x) && (x>0);
    valid_ICs = @(x) all(size(initial_guess_ICs==[6 1]));
    valid_f_handle = @(x) isa(x, 'function_handle');
    valid_bool = @(x) isa(x, 'logical');
    valid_int = @(x) mod(x,1)==0;
    
    addRequired(p,'free_vars_guess',@isnumeric);
    addRequired(p,'constraint_func',valid_f_handle);
    
    addParameter(p,'jacobian_func',default_jacobian_func,valid_f_handle);
    addParameter(p,'bool_verbose',default_bool_verbose,valid_bool);
    addParameter(p,'bool_timer',default_bool_timer,valid_bool);
    addParameter(p,'tol',default_tol,valid_scalar_pos_num);
    addParameter(p,'max_iter',default_max_iter,valid_int);
    addParameter(p,'bool_full_hist',default_bool_full_hist,valid_bool);
    
    parse(p,free_vars_guess,constraint_func,varargin{:});
    
    free_vars_guess = p.Results.free_vars_guess;
    constraint_func = p.Results.constraint_func;

    jacobian_func = p.Results.jacobian_func;
    bool_verbose = p.Results.bool_verbose;
    bool_timer = p.Results.bool_timer;
    bool_full_hist = p.Results.bool_full_hist;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    
    if ~jacobian_func(free_vars_guess)
        bool_jacobian = false;
        %jacobian_func = @numerical_jacobian;
    else
        bool_jacobian = true;
    end
 
    %% Setup

    % Build initial big X
    chi = free_vars_guess;
    
    % Initial constraint vector, chi is free variable vector
    g = constraint_func(chi);

    %% Multiple Shooting Loop
    
    % Stop condition for loop
    bool_converged = false;
    
    % List of potential scalar multiplicative factors for shooting update
    %attenuation_factors = linspace(0.1,1,5);
    attenuation_factors = (1);
    
    %err_mag = 100;
    err_mag = norm(g);
    err_mag_old = err_mag;
    err_mag_hist = [err_mag];
    
    if bool_full_hist
        chi_hist = [chi];
        if bool_jacobian
            jac_hist(:,:,1) = jacobian_func(chi);
        else
            jac_hist(:,:,1) = numerical_jacobian(constraint_func,chi);
        end
    end

    counter = 0;
    
    if bool_verbose
        fprintf("Starting shooting algorithm. ")
    end
    if bool_jacobian
        fprintf("Using supplied analytical jacobian.\n")
    else
        fprintf("Using numerical jacobian.\n")
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
        fprintf("Computing jacobian...")
        jacobian_tic = tic;
        if bool_jacobian
            dgdchi = jacobian_func(chi);
        else
            dgdchi = numerical_jacobian(constraint_func,chi);
        end
        fprintf("done.\n")
        toc(jacobian_tic)
        
        if det(dgdchi*dgdchi')==0
            error("Jacobian not invertible.")
        end

        % Test updates with different values of gamma to find greatest error reduction
        err_mag_vec = NaN(1,length(attenuation_factors));
        gtest = g;
        inv_matrix = dgdchi'/(dgdchi*dgdchi'); % precompute this for speed
        parfor c = 1:length(attenuation_factors)
            delta_chi = inv_matrix*(-attenuation_factors(c)*g);
            chi_test = chi + delta_chi;

            gtest = constraint_func(chi_test);
            err_mag = norm(gtest);
            err_mag_vec(c) = err_mag;
        end
        
        % Find correction that reduces error the most
        min_atten_factor_index = find(err_mag_vec == min(err_mag_vec));
        err_mag = err_mag_vec(min_atten_factor_index);
        
        % If error did not decrease, stop loop
%         if err_mag >= err_mag_old
%             err_mag_hist = [err_mag_hist, err_mag_old];
%             fprintf("Residual did not decrease. Terminating with a final residual of %d.\n", err_mag_old)
%             break
%         else
%             err_mag_old = err_mag; % If error continues to decrease, reset old error for comparison
%         end
        
        % Use gamma corresponding to maximum error decrease
        attenuation_factor = attenuation_factors(min_atten_factor_index);

        % Compute update to free variable vector chi
        delta_chi = inv_matrix*(-attenuation_factor*g);
        
        % Update free variables
        chi = chi + delta_chi;
        
        % Compute constraint vector
        g = constraint_func(chi);
        
        % Compute error
        err_mag = norm(g);
        
        % Save error magnitude
        err_mag_hist = [err_mag_hist, err_mag];
        % Save other history if desired
        if bool_full_hist
            chi_hist = [chi_hist,chi];
            %jac_hist = cat(3,jac_hist,dgdchi);
        end
        
        counter = counter+1;

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
    
    %% Save results
    if bool_full_hist
        results = struct('free_vars',chi,'err_hist',err_mag_hist,'iterations',counter,'chi_hist',chi_hist,'jac_hist',jac_hist);
    else
        results = struct('free_vars',chi,'err_hist',err_mag_hist,'iterations',counter);
    end
end