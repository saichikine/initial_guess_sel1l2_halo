%% Low-Thrust Transfer Construction Test
% Use multiple shooting to correct transfer with low-thrust enabled arcs
% Use variable direction and include possibility for coast arcs


%% Wipeout HD
clear; close all; clc;
set(groot,'defaultLineLineWidth', 1.5)
set(groot,'defaultAxesFontSize', 16)
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
p = gcp

%% Load stuff

load('EM_sys_stuff.mat');

eps_event = 5e-2; % epsilon bound for node placement event

G = 6.67408e-20;
T_EM = 2*pi*sqrt(L_EM^3/(moon.mu+earth.mu));

% Spacecraft params
m_sc = 200; % 200kg
D_sc = 1.6; % 1.2m diameter
Tmax = 5; % newtons
exh_vel_dim = 25000; % m/s

% Normalization stuff (divide by these to get normalized versions)
T_EM = 2*pi*sqrt(L_EM^3/(moon.mu+earth.mu));
DU = L_EM;
TU = 1/(2*pi/T_EM);
VU = DU/TU;
AU = DU/TU^2;
MU = m_sc;
FU = MU*DU/TU^2;
normalizers = struct('time_norm',TU,'vel_norm',VU,'accel_norm',AU,'force_norm',FU','m_norm',m_sc);

% Normalize quantities of interest
m_sc_normalized = m_sc/MU;
exh_vel = exh_vel_dim/VU;

%%
ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

initial_state = [0.870183560551714
  -0.059444427060739
                   0
  -0.010471778551178
  -0.175136544077527
                   0];
               
target_state = [1.115597277192111
  -0.056398494349488
                   0
  -0.008555007719687
   0.157211765557303
                   0];

%% Simulate and plot forward and backward trajectories

phase1_num_stages = 38;
phase2_num_stages = 37;

phase1_time = 2.55;
phase2_time = 3.9;
fprintf('Simulating initial guesses for departure and arrival phases...')
tic
[phase1_t_hist,phase1_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_EM), linspace(0, phase1_time, phase1_num_stages), initial_state, ode_opts);
[phase2_t_hist,phase2_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_EM), linspace(0, -phase2_time, phase2_num_stages), target_state, ode_opts);
fprintf('done.\n')
toc

% Plot
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot(x_L1, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', '$$L_1$$'); hold on % L1 location
plot(x_L2, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', '$$L_2$$'); hold on % L2 location
plot(phase1_state_hist(1,1), phase1_state_hist(1,2), 'ok', 'markerfacecolor','g','DisplayName','Departure Phase Initial Point'); hold on
plot(phase2_state_hist(1,1), phase2_state_hist(1,2), 'ok', 'markerfacecolor',[244,179,66]./255, 'DisplayName', 'Arrival Phase Initial Point'); hold on
plot(phase1_state_hist(:,1), phase1_state_hist(:,2), 'g-','DisplayName','Departure Phase Trajectory'); hold on
plot(phase2_state_hist(:,1), phase2_state_hist(:,2), 'LineStyle','-','Color',[244,179,66]./255,'DisplayName', 'Arrival Phase Trajectory'); hold on
title('Initial Guess Departure and Arrival Phases')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
legend();
axis equal
hold off

%% Define initial guess for multiple shooting

% dep_node_index = [1; 60; 120; 180; 248; 265; 282];
% arr_node_index = [1; 55; 110; 170; 232; 249; 268; size(X_hist_arr_guess,1)];
%phase1_indices = floor(linspace(1,length(phase1_state_hist)-1,phase1_num_stages));
%phase2_indices = floor(linspace(1,length(phase2_state_hist),phase2_num_stages));
phase1_indices = 1:1:phase1_num_stages-1;
phase2_indices = 1:1:phase2_num_stages;

% Column vectors going across (before assembly into column vector)
% Don't put in end of departure trajectory, just use beginning of next
% phase as the next node
phase1_stage_states = NaN(6,length(phase1_indices));
phase2_stage_states = NaN(6,length(phase2_indices));

% Assemble departure leg nodes
for i = 1:length(phase1_indices)
    phase1_stage_states(:,i) = phase1_state_hist(phase1_indices(i),:)';
end

% Assemble arrival leg nodes
for i = 1:length(phase2_indices)
    phase2_stage_states(:,i) = phase2_state_hist(phase2_indices(i),:)';
end
phase2_stage_states = fliplr(phase2_stage_states); % Flip left to right since these are propagated backwards in time

% Assemble departure and arrival nodes
stage_states = [phase1_stage_states, phase2_stage_states];

% Number of nodes
N = size(stage_states,2);

% Add on initial mass to all node states
stage_states = [stage_states; m_sc*ones(1,N)]; % initial guess is all nodes have starting mass of initial mass m0 at start of trajectory

% Thrust mvectors (ith thrust vec is thrust vec over stage i)
nom_thrust_mN = 0.5; % nominal thrust [mN]
stage_thrust_vecs = repmat([nom_thrust_mN/1000/FU;0;0],1,N-1);
stage_thrust_vecs(:,length(stage_thrust_vecs)+1) = zeros(3,1);
if length(stage_thrust_vecs) ~= N
    error("Too many control vectors specified; should have %i values, but has %i.\n",N,length(stage_thrust_vecs))
end

% Append control vectors to states
stage_states = [stage_states; stage_thrust_vecs];

% Integration times between stages
phase1_stage_times = NaN(length(phase1_indices)+1,1);
for i = 1:length(phase1_stage_times)-1
    phase1_stage_times(i) = phase1_t_hist(phase1_indices(i));
end
phase1_stage_times(end) = phase1_t_hist(end);

phase2_stage_times = NaN(length(phase2_indices),1);
for i = 1:length(phase2_stage_times)
    phase2_stage_times(i) = phase2_t_hist(phase2_indices(i));
end
phase2_stage_times = flipud(phase2_stage_times);

stage_integration_times = [abs(diff(phase1_stage_times)); abs(diff(phase2_stage_times))]./1;
                        
% Assemble free variable vector initial guess
%slack_vars_initial_guess = zeros(size(stage_integration_times));
slack_vars_initial_guess = (stage_integration_times-0.001).^(1/2);
chi_initial_guess = [reshape(stage_states,[],1); stage_integration_times; slack_vars_initial_guess];

% Create thrust/no-thrust flag vector
thrust_flags = ones(N-1,1);
% thrust_flags(1) = 0;
% thrust_flags(end) = 0;
% thrust_flags = zeros(N-1,1);
% thrust_flags(6:10) = 1;
% 
% stage_integration_times(5) = 0.25;
% stage_integration_times(6) = 0.3;

%% Plot with nodes highlighted

arc_initial_states = stage_states;
arc_integration_times = stage_integration_times;

X_hist_total_guess = [];
u_hist_total_guess = [];

fprintf("Simulating initial guess trajectory...")
tic
for i = 1:N-1
        
    if thrust_flags(i) == 0
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_EM,exh_vel), [0, arc_integration_times(i)], [arc_initial_states(1:7,i); zeros(3,1)], ode_opts);
    else
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_EM,exh_vel), [0, arc_integration_times(i)], arc_initial_states(:,i), ode_opts);
    end
    
    % Go back through and compute control history
    for j = 1:length(t_hist_arc)
        u_hist_total_guess = [u_hist_total_guess; X_hist_arc(j,8:end)];
    end
    % Save total state history
    X_hist_total_guess = [X_hist_total_guess; X_hist_arc];
end
fprintf("done.\n")
toc

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', '$$L_1$$'); hold on % L1 location
plot3(x_L2, 0, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', '$$L_2$$'); hold on % L2 location
for i = 1:size(arc_initial_states,2)
    node_string = "Node" + num2str(i);
    plot3(arc_initial_states(1,i), arc_initial_states(2,i), arc_initial_states(3,i), 'o', 'DisplayName', node_string); hold on
end
%plot(X_hist_total_guess(1,1), X_hist_total_guess(1,2), X_hist_total_guess(1,3), 'ok', 'markerfacecolor',[244,179,66]./255, 'DisplayName', 'Initial Point'); hold on
scatter3(X_hist_total_guess(:,1), X_hist_total_guess(:,2), X_hist_total_guess(:,3), 'r.','DisplayName', 'Low-Thrust Transfer Trajectory'); hold on
quiver3(X_hist_total_guess(:,1),X_hist_total_guess(:,2),X_hist_total_guess(:,3),u_hist_total_guess(:,1),u_hist_total_guess(:,2),u_hist_total_guess(:,3), 1.1,'DisplayName', 'Control Accel Vectors'); hold on
title('Initial Guess Low-Thrust Transfer')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
%legend();
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
% for i = 1:length(X_hist)
%     scatter3(X_hist_total_guess(i,1), X_hist_total_guess(i,2), X_hist_total_guess(i,3), 'r.'); hold on
%     drawnow
%     pause(0.001/norm(X_hist(i,4:6)))
% end
hold off

%% Multiple shooting (endpoint constraints)
initial_state_full = stage_states(:,1);
target_state_full = stage_states(:,end);
constraint_func = @(chi) LT_continuity_endpoint_constraint_vec(chi,N,thrust_flags,initial_state_full,target_state_full,mu_EM,exh_vel);
fprintf("Doing multiple shooting...\n\n")
tic
ms_results = shooter(chi_initial_guess,constraint_func,'tol',1e-12,'max_iter',300,'bool_verbose',true,'bool_timer',true,'bool_full_hist',true);
fprintf("done.\n\n")
toc
%% Compute Final Trajectory Results
nX = 10;
chi = ms_results.free_vars;
arc_initial_states = reshape(chi(1:10*N),10,[]);
arc_integration_times = chi(nX*N+1:nX*N+N);
slack_vars = chi(nX*N+N+1:end);

X_hist_total = [];
u_hist_total = [];
t_hist_total = [];
t_last = 0;

fprintf("Simulating converged trajectory...")
tic
for i = 1:N-1
      
    if thrust_flags(i) == 0
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_EM,exh_vel), [0+t_last, arc_integration_times(i)+t_last], [arc_initial_states(1:7,i); zeros(3,1)], ode_opts);
    else
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_EM,exh_vel), [0+t_last, arc_integration_times(i)+t_last], arc_initial_states(:,i), ode_opts);
    end
    t_last = t_hist_arc(end);
    
    % Go back through and compute control history
    for j = 1:length(t_hist_arc)
        u_hist_total = [u_hist_total; X_hist_arc(j,8:end)];
    end
    % Save total state history
    X_hist_total = [X_hist_total; X_hist_arc];
    t_hist_total = [t_hist_total; t_hist_arc];
end
fprintf("done.\n")
toc

% Plot control history
figure
addToolbarExplorationButtons(gcf)
hold on
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,1).*FU*1000, '.-', 'DisplayName', 'x control accel');
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,2).*FU*1000, '.-', 'DisplayName', 'y control accel');
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,3).*FU*1000, '.-', 'DisplayName', 'z control accel');
hold off
xlabel('Time [days]')
ylabel('Control Thrust [$$mN$$]')
legend()

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', '$$L_1$$'); hold on % L1 location
plot3(x_L2, 0, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', '$$L_2$$'); hold on % L2 location
for i = 1:size(arc_initial_states,2)
    node_string = "Node" + num2str(i);
    %plot3(arc_initial_states(1,i), arc_initial_states(2,i), arc_initial_states(3,i), 'o', 'DisplayName', node_string); hold on
end
%plot(X_hist_total_guess(1,1), X_hist_total_guess(1,2), X_hist_total_guess(1,3), 'ok', 'markerfacecolor',[244,179,66]./255, 'DisplayName', 'Initial Point'); hold on
scatter3(initial_state_full(1), initial_state_full(2), initial_state_full(3), 'x', 'DisplayName', 'Initial Position');
scatter3(target_state_full(1), target_state_full(2), target_state_full(3), 'x', 'DisplayName', 'Target Position');
scatter3(X_hist_total(:,1), X_hist_total(:,2), X_hist_total(:,3), 'r.','DisplayName', 'Low-Thrust Transfer Trajectory'); hold on
quiver3(X_hist_total(:,1),X_hist_total(:,2),X_hist_total(:,3),u_hist_total(:,1),u_hist_total(:,2),u_hist_total(:,3), 1.1,'DisplayName', 'Control Accel Vectors'); hold on
title('Final Low-Thrust Transfer')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
legend();
h = get(gca,'DataAspectRatio');
% if h(3)==1
%       set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
% else
%       set(gca,'DataAspectRatio',[1 1 h(3)])
% end
% for i = 1:length(X_hist_total)
%     scatter3(X_hist_total(i,1), X_hist_total(i,2), X_hist_total(i,3), 'r.'); hold on
%     drawnow
%     pause(0.001/norm(X_hist_total(i,4:6)))
% end
hold off
%%
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot3(arc_initial_states(1,1),arc_initial_states(2,1),arc_initial_states(3,1),'g.','DisplayName','Target Initial Point'); hold on
plot3(arc_initial_states(1,end),arc_initial_states(2,end),arc_initial_states(3,end),'r.','DisplayName','Target Final Point'); hold on
for i = 1:size(arc_initial_states,2)
    node_string = "Node" + num2str(i);
    plot3(arc_initial_states(1,i), arc_initial_states(2,i), arc_initial_states(3,i), 'o', 'DisplayName', node_string); hold on
end
plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor',[244,179,66]./255, 'DisplayName', 'Initial Point'); hold on
plot3(X_hist_total(1:end,1), X_hist_total(1:end,2), X_hist_total(1:end,3), 'LineStyle','-','Color',[244,179,66]./255,'DisplayName', 'Low-Thrust Transfer Trajectory'); hold on
%quiver3(X_hist_total(1,1),X_hist_total(1,2),X_hist_total(1,3),X_hist_total(1,4),X_hist_total(1,5),X_hist_total(1,6),2,'DisplayName', 'Velocity Vectors'); hold on
quiver3(X_hist_total(1:end,1),X_hist_total(1:end,2),X_hist_total(1:end,3),u_hist_total(1:end,1),u_hist_total(1:end,2),u_hist_total(1:end,3),'DisplayName', 'Control Accel Vectors'); hold off
title('Converged Low-Thrust Transfer')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
legend();
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
% for i = 1:length(X_hist_total)
%     scatter3(X_hist_total(i,1), X_hist_total(i,2), X_hist_total(i,3), 'r.'); hold on
%     %quiver3(X_hist_total(i,1),X_hist_total(i,2),X_hist_total(i,3),u_hist_total(i,1),u_hist_total(i,2),u_hist_total(i,3),'DisplayName', 'Control Accel Vectors'); hold on
%     drawnow
%     %pause(0.00001/norm(X_hist_total(i,4:6)))
%     pause(0.0001)
% end

%% Plot mass history and control history

figure
addToolbarExplorationButtons(gcf)
plot(t_hist_total.*time_norm./60./60./24,(200-X_hist_total(:,7)).*1000,'DisplayName','Grams of fuel used')
title('Fuel Used vs Time')
xlabel('Time [days]')
ylabel('Fuel Mass Used [g]')
grid on
legend();

figure
addToolbarExplorationButtons(gcf)
subplot(3,1,1)
plot(t_hist_total.*time_norm./60./60./24,u_hist_total(:,1).*X_hist_total(:,7).*1000,'DisplayName','x component');
title('Componentwise Control History vs Time')
ylabel('Thrust [mN]')
grid on
legend();
subplot(3,1,2)
plot(t_hist_total.*time_norm./60./60./24,u_hist_total(:,2).*X_hist_total(:,7).*1000,'DisplayName','y component');
ylabel('Thrust [mN]')
grid on
legend();
subplot(3,1,3)
plot(t_hist_total.*time_norm./60./60./24,u_hist_total(:,3).*X_hist_total(:,7).*1000,'DisplayName','z component');
xlabel('Time [days]')
ylabel('Thrust [mN]')
grid on
legend();