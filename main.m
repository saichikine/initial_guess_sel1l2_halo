% Halo connection test
% Sun-Earth system

%% Apathy is death
clear; close all; clc;
set(groot,'defaultLineLineWidth', 1.5)
set(groot,'defaultAxesFontSize', 16)
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
p = gcp;

%% character creation

earth = celestial_body;
earth.radius = 6378.13643; % [km]
earth.mu = 3.986e5; % [km^3/s^2]
 
moon = celestial_body;
moon.radius = 1738; % [km]
moon.mu = 4902.799; % [km^3/s^2]

sun = celestial_body;
sun.radius = 696000; % [km]
sun.mu = 1.327e11; % [km^3/s^2]

mu_ME = moon.mu/(moon.mu + earth.mu);
mu_SE = earth.mu/(earth.mu + sun.mu);

L_EM = 384400; % [km], Earth-Moon distance
L_SE = 149598023; % [km], Sun-Earth distance

L_points = lagrangePoints(mu_SE);
x_L1 = L_points(1,1);
x_L2 = L_points(1,2);

G = 6.67408e-20;

%% Define spacecraft
g0 = 9.80665; % [m/s^2] g0 for Isp

% 6U cubesat (from ians paper)
cubesat = spacecraft;
spacecraft.mass = 14; % [kg]
spacecraft.max_thrust = 0.4; % [mN]
spacecraft.Isp = 1250; % [s]
spacecraft.Ve = spacecraft.Isp*g0; % [m/s]

%%
% Normalization stuff (divide by these to get normalized versions)
LU = L_SE;
T_SU = 2*pi*sqrt(LU^3/(earth.mu+sun.mu));
DU = L_SE; % [km]
TU = 1/(2*pi/T_SU); % [s]
VU = DU/TU; % [km/s]
AU = DU/TU^2; % [km/s^2]
MU = spacecraft.mass; % [kg]
FU = MU*DU/TU^2; % [kg*km/s^2] = [kN]
normalizers = struct('time_norm',TU,'vel_norm',VU,'accel_norm',AU,'force_norm',FU','m_norm',MU);

% Normalize quantities of interest
m_sc = spacecraft.mass/MU;
exh_vel = spacecraft.Ve/1000/VU; % convert to km/s, then to nondimensional velocities
max_thrust_mag = spacecraft.max_thrust/1000/1000/FU; % convert to kN, then to nondimensional force units

% Find ICs for L1 and L2 Lyapunov orbits with same Jacobi Constant
% Use bisection to find orbit with matching Jacobi constant

% Numerical parameters
eps = 5e-2; % epsilon bound for node placement event

%% First (reference) halo orbit
L_km = L_SE;
Az_L1 = 200000; %[km]

halo_IC_L1_southern = halo_computeplot(mu_SE, Az_L1, L_km, "L1", "south", 0);
fprintf("First orbit found.\n")

halo1_JC = jacobi_constant(halo_IC_L1_southern{1},mu_SE);
fprintf("First orbit's Jacobi constant is %d.\n", halo1_JC);
%% Find second orbit to match Jacobi constant with first orbit

% initial interval
a = Az_L1/2;
b = Az_L1*2;

% numerical setup
bi_count = 0;
max_bi_count = 30;
bi_tol = 1e-8; % lower for SE system

zero_func = @(Az) zero_func_params(Az, mu_SE, L_km, halo1_JC);
tic
while bi_count < max_bi_count

    c = (a+b)/2;
    
    if abs(zero_func(c)) < bi_tol || (b-a)/2 < bi_tol/100
        sln = c;
        break;
    end
    
    bi_count = bi_count+1;
    fprintf('Iteration count: %i\n',bi_count);
    
    if sign(zero_func(c)) == sign(zero_func(a))
        a = c;
    else
        b = c;
    end
end
toc
fprintf("Done, with %i iterations.\n",bi_count);
%% Show Jacobi constant results
clc;

Az_L2 = c;
halo_IC_L2_southern = halo_computeplot(mu_SE, Az_L2, L_km, "L2", "north", 0);
halo2_JC = jacobi_constant(halo_IC_L2_southern{1},mu_SE);
fprintf('Found Jacobi constant: %d\n',halo2_JC);
fprintf('Target Jacobi constant: %d\n',halo1_JC);

%% Integrate and plot both orbits

ode_opts = odeset('RelTol',5e-14,'AbsTol',1e-15);

T_L1 = halo_IC_L1_southern{2};
X0_L1 = halo_IC_L1_southern{1};
[~, X_hist_L1] = ode113(@(t,X) CR3BP(t,X,mu_SE), [0 T_L1], X0_L1, ode_opts);

T_L2 = halo_IC_L2_southern{2};
X0_L2 = halo_IC_L2_southern{1};
[~, X_hist_L2] = ode113(@(t,X) CR3BP(t,X,mu_SE), [0 T_L2], X0_L2, ode_opts);

figure
addToolbarExplorationButtons(gcf)
hold on
scatter3(x_L1*L_km, 0, 0, 'd', 'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','L1'); hold on
scatter3(x_L2*L_km, 0, 0, 'd', 'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','L2');
plot3((1-mu_SE)*L_km, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 10, 'DisplayName', 'Earth'); hold on % Smaller primary
plot3(X_hist_L1(:,1)*L_km, X_hist_L1(:,2)*L_km, X_hist_L1(:,3)*L_km, 'b-','DisplayName', '$$L_1$$ Orbit'); hold on
plot3(X_hist_L2(:,1)*L_km, X_hist_L2(:,2)*L_km, X_hist_L2(:,3)*L_km, 'r-','DisplayName', '$$L_2$$ Halo Orbit'); hold on
hold off
grid on
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1])
else
      set(gca,'DataAspectRatio',[1 1 1])
end
legend()
xlabel('$$x [km]$$')
ylabel('$$y [km]$$')
zlabel('$$z [km]$$')
view([-37.1298828125,29.7275390625]);
set(gcf,'color','w')
%title("Initial and target halo orbits")

%% Poincare maps for both orbits
clc
fprintf('Computing Poincare maps...\n')
tic
map1 = poincare_map(mu_SE,L_km,TU,halo_IC_L1_southern,200*24*60*60,"unstable","interior",1,100);
map2 = poincare_map(mu_SE,L_km,TU,halo_IC_L2_southern,200*24*60*60,"stable","exterior",2,100);
toc
fprintf('done.\n')
%% Plot Poincare maps

figure
addToolbarExplorationButtons(gcf);
plot3(map1.map_points(1,:),map1.map_points(2,:),map1.map_points(3,:),'r.-','DisplayName','Unstable manifold from L1 orbit'); hold on
plot3(map2.map_points(1,:),map2.map_points(2,:),map2.map_points(3,:),'b.-','DisplayName','Stable manifold to L2 orbit'); hold off
% xlim([-0.05 0.05])
% zlim([0 0.3])
view([90 0])
% xlabel('$$\dot{x}$$')
% ylabel('$$\dot{y}$$')
% zlabel('$$\dot{z}$$')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
title('Combined Poincare map')
legend()
set(gcf,'color','w')
 
%%
% Now use multiple shooting to connect the two manifold segments

% Pick trajectories that don't quite intersect on Poincare map but are
% close
crit_index1 = 7; % manifold trajectory index based on looking at Poincare maps
crit_index2 = 56;

% crit_index1 = 198;
% crit_index2 = 97;

initial_state = map1.man_trajs{crit_index1}{2}(1:6,1);
target_state = map2.man_trajs{crit_index2}{2}(1:6,1);

%% Simulate forward and backward trajs

ode_opts = odeset('RelTol',5e-14,'AbsTol',1e-20);

phase1_time = 2.15;
phase2_time = 2.12;

% phase1_time = 2.1;
% phase2_time = 2.15;
fprintf('Simulating initial guesses for departure and arrival phases...')
tic
% [phase1_t_hist,phase1_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_SE), linspace(0, phase1_time, phase1_num_stages), initial_state, ode_opts);
% [phase2_t_hist,phase2_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_SE), linspace(0, -phase2_time, phase2_num_stages), target_state, ode_opts);
[phase1_t_hist,phase1_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_SE), [0 phase1_time], initial_state, ode_opts);
[phase2_t_hist,phase2_state_hist] = ode113(@(t,X) CR3BP(t,X,mu_SE), [0 -phase2_time], target_state, ode_opts);
fprintf('done.\n')
toc

% Plot
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
scatter3(x_L1, 0, 0, 'd', 'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','L1'); hold on
scatter3(x_L2, 0, 0, 'd', 'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','L2');
plot3(1-mu_SE, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 10, 'DisplayName', 'Earth'); hold on % Smaller primary
plot3(phase1_state_hist(1,1), phase1_state_hist(1,2), phase1_state_hist(1,3), 'ok', 'markerfacecolor','g','DisplayName','Departure Phase Initial Point'); hold on
plot3(phase2_state_hist(1,1), phase2_state_hist(1,2), phase2_state_hist(1,3), 'ok', 'markerfacecolor',[244,179,66]./255, 'DisplayName', 'Arrival Phase Initial Point'); hold on
plot3(phase1_state_hist(:,1), phase1_state_hist(:,2), phase1_state_hist(:,3), 'g-','DisplayName','Departure Phase Trajectory'); hold on
plot3(phase2_state_hist(:,1), phase2_state_hist(:,2), phase2_state_hist(:,3), 'LineStyle','-','Color',[244,179,66]./255,'DisplayName', 'Arrival Phase Trajectory'); hold on
title('Initial Guess Departure and Arrival Phases')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
legend();
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1])
else
      set(gca,'DataAspectRatio',[1 1 1])
end
hold off

%% Define initial guess for multiple shooting

phase1_num_stages = 50;
phase2_num_stages = 50;
phase1_indices = floor(linspace(1,length(phase1_state_hist)-1,phase1_num_stages));
phase2_indices = floor(linspace(1,length(phase2_state_hist),phase2_num_stages));
% phase1_indices = 1:1:phase1_num_stages-1;
% phase2_indices = 1:1:phase2_num_stages;

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
nom_thrust_frac = 1e-1; % nominal thrust (between 0 and 1, multiplies Tmax)
stage_control_vecs = NaN(3,N);
for i = 1:N-1
    stage_control_vecs(:,i) = nom_thrust_frac*stage_states(4:6,i)/norm(stage_states(4:6,i));
end
stage_control_vecs(:,end) = zeros(3,1);
% stage_thrust_vecs = repmat([nom_thrust_mN/1000/FU;0;0],1,N-1);
% stage_thrust_vecs(:,length(stage_thrust_vecs)+1) = zeros(3,1);
if length(stage_control_vecs) ~= N
    error("Too many control vectors specified; should have %i values, but has %i.\n",N,length(stage_control_vecs))
end

% Append control vectors to states
stage_states = [stage_states; stage_control_vecs];

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
slack_vars_initial_guess = (stage_integration_times-5e-2).^(1/2);
chi_initial_guess = [reshape(stage_states,[],1); stage_integration_times; slack_vars_initial_guess];

% Create thrust/no-thrust flag vector
thrust_flags = ones(N-1,1);
thrust_flags(1) = 0;
thrust_flags(end) = 0;
% thrust_flags = zeros(N-1,1);
% thrust_flags(8:12) = ones(5,1);

initial_state_full = stage_states(:,1);
target_state_full = stage_states(:,end);

%% Plot with nodes highlighted

nX = 10;
chi = chi_initial_guess;
% arc_initial_states = stage_states;
% arc_integration_times = stage_integration_times;
arc_initial_states = reshape(chi(1:10*N),10,[]);
arc_integration_times = chi(nX*N+1:nX*N+N-1);
slack_vars = chi(nX*N+N:end);

X_hist_total_guess = [];
u_hist_total_guess = [];

fprintf("Simulating initial guess trajectory...")
tic
for i = 1:N-1
        
    if thrust_flags(i) == 0
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_SE,exh_vel,max_thrust_mag), [0, arc_integration_times(i)], [arc_initial_states(1:7,i); zeros(3,1)], ode_opts);
    else
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_SE,exh_vel,max_thrust_mag), [0, arc_integration_times(i)], arc_initial_states(:,i), ode_opts);
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
plot3(1-mu_SE, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
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
      set(gca,'DataAspectRatio',[1 1 1])
else
      set(gca,'DataAspectRatio',[1 1 1])
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
constraint_func = @(chi) LT_continuity_endpoint_constraint_vec(chi,N,thrust_flags,initial_state_full,target_state_full,mu_SE,exh_vel,max_thrust_mag);
fprintf("Doing multiple shooting...\n\n")
tic
ms_results = shooter(chi_initial_guess,constraint_func,'tol',1e-14,'max_iter',300,'bool_verbose',true,'bool_timer',true,'bool_full_hist',false);
fprintf("done.\n\n")
toc
%% Compute Final Trajectory Results
nX = 10;
chi = ms_results.free_vars;
arc_initial_states = reshape(chi(1:10*N),10,[]);
arc_integration_times = chi(nX*N+1:nX*N+N-1);
slack_vars = chi(nX*N+N:end);

ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-19);

X_hist_total = [];
u_hist_total = [];
t_hist_total = [];
t_last = 0;

fprintf("Simulating converged trajectory...")
tic
for i = 1:N-1
      
    if thrust_flags(i) == 0
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_SE,exh_vel,max_thrust_mag), [0+t_last, arc_integration_times(i)+t_last], [arc_initial_states(1:7,i); zeros(3,1)], ode_opts);
    else
        [t_hist_arc,X_hist_arc] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_SE,exh_vel,max_thrust_mag), [0+t_last, arc_integration_times(i)+t_last], arc_initial_states(:,i), ode_opts);
    end
    t_last = t_hist_arc(end);
    
    % Go back through and compute control history
    for j = 1:length(t_hist_arc)
        u_hist_total = [u_hist_total; max_thrust_mag*FU*1000*1000*X_hist_arc(j,8:end)];
    end
    % Save total state history
    X_hist_total = [X_hist_total; X_hist_arc];
    t_hist_total = [t_hist_total; t_hist_arc];
end
[t_hist_targ_orb, X_hist_targ_orb] = ode113(@(t,X) CR3BP_cart_control(t,X,mu_SE,exh_vel,max_thrust_mag), [0+t_last, T_L2+t_last], [arc_initial_states(1:7,end); zeros(3,1)], ode_opts);
for i = 1:length(t_hist_targ_orb)
    u_hist_total = [u_hist_total; zeros(3,1)'];
end
X_hist_total = [X_hist_total; X_hist_targ_orb];
t_hist_total = [t_hist_total; t_hist_targ_orb];
fprintf("done.\n")
toc

% Plot control history
figure
addToolbarExplorationButtons(gcf)
hold on
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,1), '.-', 'DisplayName', 'x thrust');
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,2), '.-', 'DisplayName', 'y thrust');
plot(t_hist_total.*TU/60/60/24, u_hist_total(:,3), '.-', 'DisplayName', 'z thrust');
hold off
xlabel('Time [days]')
ylabel('Control Thrust [$$mN$$]')
legend()
grid on

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_SE, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 10, 'DisplayName', 'Earth'); hold on % Smaller primary
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
quiver3(X_hist_total(:,1),X_hist_total(:,2),X_hist_total(:,3),u_hist_total(:,1),u_hist_total(:,2),u_hist_total(:,3), 1.1,'DisplayName', 'Thrust Vectors'); hold on
title('Final Low-Thrust Transfer')
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
%legend();
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1])
else
      set(gca,'DataAspectRatio',[1 1 1])
end
% for i = 1:length(X_hist_total)
%     scatter3(X_hist_total(i,1), X_hist_total(i,2), X_hist_total(i,3), 'r.'); hold on
%     drawnow
%     pause(0.001/norm(X_hist_total(i,4:6)))
% end
hold off

%% Save results to file

ms_results_initial_guess = struct('ms_results',ms_results,'chi',chi,'arc_initial_states',arc_initial_states,'arc_integration_times',arc_integration_times,'slack_vars',slack_vars,'num_nodes',N,'thrust_flags',thrust_flags,'L1_start_orbit',halo_IC_L1_southern,'L2_end_orbit',halo_IC_L2_southern);
save('halo_trans_ig_100_good.mat','ms_results_initial_guess');
%%
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_SE, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 10, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
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
      set(gca,'DataAspectRatio',[1 1 1])
else
      set(gca,'DataAspectRatio',[1 1 1])
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

%% Inline functions

% function we are trying to make zero
function zero_val = zero_func_params(Az, mu_SE, L, jacobi_target)
    halo_IC = halo_computeplot(mu_SE, Az, L, "L2", "north", 0);
    halo_JC = jacobi_constant(halo_IC{1},mu_SE);
    zero_val = halo_JC - jacobi_target;
end