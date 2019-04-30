% lyapunov connection test

%% Apathy is death
clear; close all; clc;

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

mu_EM = moon.mu/(moon.mu + earth.mu);
mu_SE = earth.mu/(earth.mu + sun.mu);

L_EM = 384400; % [km], Earth-Moon distance
L_SE = 149598023; % [km], Sun-Earth distance

L_points = lagrangePoints(mu_EM);
x_L1 = L_points(1,1);
x_L2 = L_points(1,2);

EM_period = 2*pi*sqrt(L_EM^3/(moon.mu+earth.mu));

prop_time = (10*24*60*60)*2*pi/EM_period;

% Find ICs for L1 and L2 Lyapunov orbits with same Jacobi Constant
% Use bisection to find orbit with matching Jacobi constant

% Numerical parameters
eps = 5e-2; % epsilon bound for node placement event

%% First (reference) orbit
L_km = L_EM;
Ax_km_L1 = 70000;
x0_L1 = Ax_km_L1/L_km;

% Set up inputs for first orbit
lin_state_L1 = @(t) linear_planar_L1L2_dynamics(t,mu_EM,1,x0_L1); % analytical sln handle
node_event_L1 = @(t,X) event_nodes_analytical(t,X,lin_state_L1,eps); % patch point event handle (places nodes/patch points when this is triggered)
L1_orbit_period = linear_planar_L1L2_orbit_period(mu_EM,1); % approximate orbit period
L1_ICs = lin_state_L1(0);

% Multiple shooting to find first (reference) orbit
L1_orbit = multiple_shooting_analytical(L1_ICs,L1_orbit_period,mu_EM,lin_state_L1,node_event_L1,'plot_flag',1);
fprintf('First orbit Jacobi constant = %d\n',L1_orbit.jacobi_constant);

%% Find second orbit to match Jacobi constant with first orbit

% initial interval
a = Ax_km_L1/2;
b = Ax_km_L1*2;

% numerical setup
bi_count = 0;
max_bi_count = 20;
bi_tol = 1e-3;

zero_func = @(Ax_km) zero_func_params(Ax_km, mu_EM, L_km, eps, L1_orbit);
tic
while bi_count < max_bi_count

    c = (a+b)/2;
    
    if abs(zero_func(c)) < bi_tol || (b-a)/2 < bi_tol
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

lin_state_L2 = @(t) linear_planar_L1L2_dynamics(t, mu_EM, 2, sln*1.01/L_km);
node_event_L2 = @(t,X) event_nodes_analytical(t,X,lin_state_L2,eps);
L2_orbit_period = linear_planar_L1L2_orbit_period(mu_EM,2);
L2_ICs = lin_state_L2(0);
fprintf('Starting multiple shooting for found matching orbit...')
tic
L2_orbit = multiple_shooting_analytical(L2_ICs,L2_orbit_period,mu_EM,lin_state_L2,node_event_L2,'plot_flag',0,'verbose_flag',0);
fprintf('done.\n')
toc
fprintf('Found Jacobi constant: %d\n',L2_orbit.jacobi_constant);
fprintf('Target Jacobi constant: %d\n',L1_orbit.jacobi_constant);

%% Poincare maps for both orbits
clc

num_man_points = 300;
fprintf('Computing Poincare maps...\n')
tic
map1 = poincare_map(mu_EM,L_km,EM_period,{L1_orbit.X0,L1_orbit.period},30*24*60*60,"unstable","exterior",1,num_man_points);
map2 = poincare_map(mu_EM,L_km,EM_period,{L2_orbit.X0,L2_orbit.period},30*24*60*60,"stable","interior",2,num_man_points);
toc
fprintf('done.\n')
%% Plot Poincare maps

figure
addToolbarExplorationButtons(gcf);
plot(map1.map_points(2,:),map1.map_points(1,:),'r.','DisplayName','Unstable manifold from L1 orbit'); hold on
plot(map2.map_points(2,:),map2.map_points(1,:),'g.','DisplayName','Stable manifold to L2 orbit'); hold off
ylim([-inf 0])
xlabel('y dot')
ylabel('y')
grid on
title('Combined Poincare map')
legend()

%% Plot initial guess
% With 300 points
crit_index1 = 131;
crit_index2 = 194;

% With 350 points
% crit_index1 = 153;
% crit_index2 = 227;

ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
X0_ig = map1.man_trajs{crit_index1}{2}(1:6,1);
X0_ig(3) = 0;
X0_ig(6) = 0;
[t_ig, X_ig] = ode113(@(t,X) CR3BP(t,X,mu_EM), [0,8], X0_ig ,ode_opts);

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot3(X_ig(1,1), X_ig(1,2), X_ig(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
scatter3(X_ig(:,1), X_ig(:,2), X_ig(:,3), 'r.','DisplayName', 'Converged Connection'); hold on
title('Multiple Shooting Converged')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;
legend();
hold off

%% Look at first orbit
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot3(map1.orbit_traj{2}(1,1), map1.orbit_traj{2}(2,1), map1.orbit_traj{2}(3,1), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
scatter3(map1.orbit_traj{2}(1,1:150), map1.orbit_traj{2}(2,1:150), map1.orbit_traj{2}(3,1:150), 'r.','DisplayName', 'First orbit'); hold on
title('Multiple Shooting Converged')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;
legend();
hold off

%% Look at second orbit

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot3(map2.orbit_traj{2}(1,1), map2.orbit_traj{2}(2,1), map2.orbit_traj{2}(3,1), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
scatter3(map2.orbit_traj{2}(1,1:150), map2.orbit_traj{2}(2,1:150), map2.orbit_traj{2}(3,1:150), 'r.','DisplayName', 'Second orbit'); hold on
title('Multiple Shooting Converged')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;
legend();
hold off
%%
% Now use multiple shooting to connect the two manifold segments
%% Build first set of nodes (departure trajectory segment)
num_orbit_nodes = 10;
node_indices1 = round(linspace(1,num_man_points,num_orbit_nodes)); % use linspace instead of colon so endpoints are included
node_points1 = [];

% Nodes along periodic orbit (with copies)
copy_num = 1;
node_points1 = [node_points1, repmat(map1.orbit_traj{2}(1:6,node_indices1(1:end-1)),1,copy_num)];

% Nodes along periodic orbit UNTIL departure point
node_points1 = [node_points1, map1.orbit_traj{2}(1:6,node_indices1(node_indices1<crit_index1))];

% Add nodes along manifold trajectory
man_nodes_count = 15;
man_index_spacing1 = floor(length(map1.man_trajs{crit_index1}{1})/man_nodes_count);
%node_points1 = [node_points1, map1.man_trajs{crit_index1}{2}(1:6,round(linspace(1,length(map1.man_trajs{crit_index1}{1}),man_nodes_count)))];
node_points1 = [node_points1, map1.man_trajs{crit_index1}{2}(1:6,1:man_index_spacing1:end)];

fprintf("Departure segment: nodes set, %i departure nodes.\n",length(node_points1))

%% Departure segment integration times

times1 = [];
times1 = [times1; repmat(diff(map1.orbit_traj{1}(node_indices1)),copy_num,1)];
times1 = [times1; diff(map1.orbit_traj{1}(node_indices1(node_indices1<crit_index1)))];
times1 = [times1; map1.man_trajs{crit_index1}{1}(1) - map1.orbit_traj{1}(node_indices1(find(node_indices1<crit_index1,1,'last')))]; % time from last orbit point to manifold traj start
% times1 = [times1; diff(map1.man_trajs{crit_index1}{1}(round(linspace(1,length(map1.man_trajs{crit_index1}{1}),man_nodes_count))))]; % manifold traj times
times1 = [times1; diff([map1.man_trajs{crit_index1}{1}(1:man_index_spacing1:end); map1.man_trajs{crit_index1}{1}(end)])]; % manifold traj times
if length(times1) ~= size(node_points1,2)
    error("Time vector is the wrong length. Should be %i, but it's %i.\n",size(node_points1,2),length(times1))
end

fprintf("Departure segment: integration times set, %i departure integration times.\n",length(times1))
%% Build second second of nodes (arrival trajectory segment)
node_indices2 = round(linspace(1,num_man_points,num_orbit_nodes));
node_points2 = [];

man_index_spacing2 = floor(length(map2.man_trajs{crit_index2}{1})/man_nodes_count);
node_points2 = fliplr([node_points2, map2.man_trajs{crit_index2}{2}(1:6,round(linspace(1,length(map2.man_trajs{crit_index2}{1}),man_nodes_count)))]);

% nodes from manifold insert point to end of that orbit
node_points2 = [node_points2, map2.orbit_traj{2}(1:6,node_indices2(node_indices2>crit_index2 & node_indices2<node_indices2(end)))]; ...
    % Don't put in last orbit state (taken care of in start of next orbit)

% Nodes along periodic orbit (3 copies, + last orbit point)
node_points2 = [node_points2, repmat(map2.orbit_traj{2}(1:6,node_indices2(1:end-1)),1,copy_num), map2.orbit_traj{2}(1:6,end)];

fprintf("Arrival segment: nodes set, %i arrival nodes.\n",length(node_points2))

%% Arrival segment integration times

times2 = [];
% no need to flip these since times are already negative)
%times2 = [times2; abs(flipud(diff(map2.man_trajs{crit_index2}{1}(round(linspace(1,length(map2.man_trajs{crit_index2}{1}),man_nodes_count))))))]; % manifold traj times
times2 = [times2; diff(flipud(map2.man_trajs{crit_index2}{1}(round(linspace(1,length(map2.man_trajs{crit_index2}{1}),man_nodes_count)))))]; % manifold traj times
times2 = [times2; map2.orbit_traj{1}(node_indices2(find(node_indices2>crit_index2,1,'first')))-map2.man_trajs{crit_index2}{1}(1)]; % time from manifold arrival point to next orbit point
times2 = [times2; diff(map2.orbit_traj{1}(node_indices2(node_indices2>crit_index2)))]; % times till end of this orbit
times2 = [times2; repmat(diff(map2.orbit_traj{1}(node_indices2)),copy_num,1)];

if length(times2) ~= size(node_points2,2)-1
    error("Time vector is the wrong length. Should be %i, but it's %i.\n",size(node_points2,2)-1,length(times2))
end

fprintf("Arrival segment: integration times set, %i arrival integration times.\n",length(times2))
    

%% Assemble free var vector (combining both segments)

%chi = [node_points1; flipud(node_points2); integration_times1; 0; flipud(integration_times2)]; % 0 in middle is integration time for connecting middle nodes
node_points = [node_points1, (node_points2)];
% times = [times1; 0; times2];
times = [times1; times2];

if length(times) ~= size(node_points,2)-1
    error("Time vector is the wrong length. Should be %i, but it's %i.\n",size(node_points,2)-1,length(times))
end

fprintf("Total nodes: %i\n",size(node_points,2));
fprintf("Total integration times: %i\n",length(times));

% set nonzero (very small) z and zdot components to zero
node_points(3,:) = zeros(1,size(node_points,2));
node_points(6,:) = zeros(1,size(node_points,2)); 

% Plot nodes
figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
scatter3(node_points1(1,:), node_points1(2,:), node_points1(3,:), 'ro','DisplayName', 'Departure Segment Nodes'); hold on
scatter3(node_points2(1,:), node_points2(2,:), node_points2(3,:), 'go','DisplayName', 'Arrival Segment Nodes'); hold on
title('Multiple Shooting Nodes')
xlabel('X')
ylabel('Y')

zlabel('Z')
grid on;
legend();
hold off

%% Multiple shooting
clc
fprintf('Doing multiple shooting...\n')
tic
connection_test = multiple_shooting_hetero(node_points, times, mu_EM, 'max_iter',50,'tol',5e-14,'plot_flag', 1, 'verbose_flag', 1);
toc
fprintf('done.\n')

%% Integrate found ICs

ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
[t_con, X_con] = ode113(@(t,X) CR3BP(t,X,mu_EM), [0 connection_test.total_time-((copy_num-1)*L2_orbit.period)], connection_test.X0, ode_opts);

%% Plot result

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot3(X_con(1,1), X_con(1,2), X_con(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
scatter3(X_con(:,1), X_con(:,2), X_con(:,3), 'r.','DisplayName', 'Converged Connection'); hold on
title('Multiple Shooting Converged')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;
legend();
hold off

%% Inline functions

% function we are trying to make zero
function zero_val = zero_func_params(Ax_km, mu_EM, L_km, eps, L1_orbit)
    lin_state_L2 = @(t) linear_planar_L1L2_dynamics(t,mu_EM,2,Ax_km/L_km);
    node_event_L2 = @(t,X) event_nodes_analytical(t,X,lin_state_L2,eps);
    L2_orbit_period = linear_planar_L1L2_orbit_period(mu_EM,2);
    L2_ICs = lin_state_L2(0);
    
    L2_orbit = multiple_shooting_analytical(L2_ICs,L2_orbit_period,mu_EM,lin_state_L2,node_event_L2,'plot_flag',0,'verbose_flag',0);
    zero_val = L2_orbit.jacobi_constant - L1_orbit.jacobi_constant;
end