function constraint_vector = LT_continuity_endpoint_constraint_vec(chi, N, thrust_flags, initial_state, final_state, mu, exh_vel, Tmax)
    
    % Computes constraint vector with variable thrust, variable direction
    % arcs
    % Accounts for coast arcs as well with thrust_flags vector
    
    %% Reformat free variable vector chi
    nX = 10;
    arc_initial_states = reshape(chi(1:nX*N),nX,[]);
    arc_integration_times = (chi(nX*N+1:nX*N+N-1));
    arc_times_slack_vars = chi(nX*N+N:end);
        
    %% Constraint vector
    g = [];
    g = [g; arc_initial_states(1:7,1)-reshape(initial_state(1:7),[],1); arc_initial_states(1:6,end)-reshape(final_state(1:6),[],1)];
    
    %% Slack variables for negative integration times
    dt_min = 5e-2;
    for i = 1:length(arc_times_slack_vars)
        g = [g; arc_integration_times(i) - arc_times_slack_vars(i)^2 - dt_min];
    end
%     negative_times_indices = find(arc_integration_times<dt_min);
%     if ~isempty(negative_times_indices)
%         for i = 1:length(negative_times_indices)
%             g = [g; arc_integration_times(negative_times_indices(i)) - dt_min];
%             %arc_integration_times(negative_times_indices(i)) = 0; 
%         end
%     end
    %% Integrate each arc forward
    
    ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
    arc_final_states = NaN(7,size(arc_initial_states,2)); 
    
    parfor i = 1:N-1
        
        if arc_integration_times(i) == 0
            arc_final_states(:,i) = arc_initial_states(1:7,i);
        else
            if thrust_flags(i) == 0
                [~,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu,exh_vel,Tmax), ([0, arc_integration_times(i)]), [arc_initial_states(1:7,i); zeros(3,1)], ode_opts);
            else
                [~,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu,exh_vel,Tmax), ([0, arc_integration_times(i)]), arc_initial_states(:,i), ode_opts);
            end

            arc_final_states(:,i) = X_hist_arc(end,1:7)';
        end
        
        % Always enforce position + velocity constraints
        g = [g; arc_final_states(:,i) - arc_initial_states(1:7,i+1)];
    end
    
    constraint_vector = g;
end