function results = poincare_map(mu, L_km, TU, orbit, prop_time_s, stability_string, int_ext_string, L_point_num, varargin)

    % Compute and plot manifold for periodic orbit
    
    %% input handling
    
    if nargin > 9
        error("Maximum of 9 arguments.")
    end
    
    if ~isstring(stability_string)
        error("Argument 6 must be a string.")
    end
    
    if ~isstring(int_ext_string)
        error("Argument 7 must be a string.")
    end
    
    stab_flag = 0;
    unstab_flag = 0;
    if strcmpi(stability_string,"stable")
        stab_flag = 1;
    elseif strcmpi(stability_string,"unstable")
        unstab_flag = 1;
    else
        error("Well this is bad.")
    end
    
    int_flag = 0;
    ext_flag = 0;
    if strcmpi(int_ext_string,"interior")
        int_flag = 1;
    elseif strcmpi(int_ext_string,"exterior")
        ext_flag = 1;
    end
    
    num_points = 100;
    if length(varargin)>0
        num_points = varargin{1};
    end

    %% it's a trap

    L_points = lagrangePoints(mu);
    xL1 = L_points(1,1);
    xL2 = L_points(1,2);
    
    prop_time = prop_time_s/TU;

    %% 

    ode_opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-20);
    X0 = [orbit{1}; reshape(eye(6),[],1)];
    [t_orbit, periodic_orbit] = ode113(@(t,X) CR3BP(t,X,mu), linspace(0, orbit{2}, num_points), X0, ode_opts);
    
    %% Manifold Stuff
    
    plane = "yz"; %temp
    dist = 1-mu;
    
    ode_opts_poincare_event = odeset('RelTol',1e-13,'AbsTol',1e-15,'Events',@(t,X) poincare_event(t,X,plane,dist,L_point_num));
    
    poincare_map_points = [];
    
    monodromy_mat = reshape(periodic_orbit(end,7:end),6,6); % STM for one orbit period

    [eig_vecs, eig_vals] = eig(monodromy_mat);

    eig_vals = diag(eig_vals);
    real_eig_vals = (eig_vals(find(imag(eig_vals)==0)));

    vS = eig_vecs(:,find(real_eig_vals==min(real_eig_vals)));
    vU = eig_vecs(:,find(real_eig_vals==max(real_eig_vals)));
    
    if stab_flag
        v = vS;
        prop_time = -prop_time;
    elseif unstab_flag
        v = vU;
    end

    figure; hold on
    addToolbarExplorationButtons(gcf)
    sp = plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
    ip = plot3(periodic_orbit(1,1), periodic_orbit(1,2), periodic_orbit(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
    l1 = plot3(xL1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
    l2 = plot3(xL2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
    po = plot3(periodic_orbit(:,1), periodic_orbit(:,2), periodic_orbit(:,3), 'b-'); hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on;
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    h = get(gca,'DataAspectRatio');
    if h(3)==1
          set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
    else
          set(gca,'DataAspectRatio',[1 1 h(3)]);
    end
    
    % Arrays to save results
    perturbed_manifold_points = NaN(6,length(t_orbit));
    poincare_map_points = NaN(6,length(t_orbit));
    % Cell array for manifold trajectories
    % Each cell element contains a time vector{1} and a state history{2}
    manifold_trajectories = cell(length(t_orbit));
    
    % Propagate points along orbit
    for i = 1:length(t_orbit)

        %eps_pos = 100/L_km;
        eps_pos = 1e-4;
        eps_vel = eps_pos/norm(periodic_orbit(i,1:3));

        eps_vec = [eps_pos*ones(3,1); eps_vel*ones(3,1)];
        
        v_t = reshape(periodic_orbit(i,7:end),6,6)*v;
        v_t = eps_vec.*v_t/norm(v_t);

        X_int = periodic_orbit(i,1:6)' + v_t;
        X_ext = periodic_orbit(i,1:6)' - v_t;
        
        if int_flag
            X0 = X_int;
        elseif ext_flag
            X0 = X_ext;
        else
            error("bad.")
        end
        
        % Enforce zero z and zdot for lyapunov case
%         X0(3) = 0;
%         X0(6) = 0; 
        
        % Save perturbed points along orbit
        perturbed_manifold_points(:,i) = X0;
        
        % propagate manifold(s) to Poincare surface of section
        [t_man, X_man, te, Xe, ie] = ode113(@(t,X) CR3BP(t,X,mu), linspace(t_orbit(i),t_orbit(i)+prop_time,100), X0, ode_opts_poincare_event);
        
        if ~isempty(Xe)
            % Save manifold trajectories
            manifold_trajectories{i}{1} = t_man;
            manifold_trajectories{i}{2} = X_man'; % transpose so states are column vecs
            poincare_map_points(:,i) = reshape(Xe,[],1);
            
%             if abs(Xe(3))>1e-24 || abs(Xe(6))>1e-12
%                 warning("Warning: nonzero z components:\n X=%d\n Xdot=%d\n",Xe(3),Xe(6));
%             end
        end
        
        % Plot manifolds
        ma = plot3(X_man(:,1), X_man(:,2), X_man(:,3), 'r-'); % interior manifold

    end
    hold off
    legend([sp, ip, l1, l2, po, ma], "Smaller Primary", "Initial Periodic Orbit Point", "L1 Point", "L2 Point", "Periodic Orbit", "Manifold");
    
    title(stability_string + " " + int_ext_string + " manifold")
    
    figure
    title('Poincare Map')
    plot(poincare_map_points(2,:), poincare_map_points(1,:), 'k.')
    xlabel('y dot')
    ylabel('y')
    grid on
    
    %% Save results
    orbit_trajectory = {t_orbit, periodic_orbit'};
    
    % all states are column vectors
    results = struct('map_points',poincare_map_points,'manifold_points',perturbed_manifold_points,...
        'man_trajs',{manifold_trajectories},'orbit_traj',{orbit_trajectory});
end