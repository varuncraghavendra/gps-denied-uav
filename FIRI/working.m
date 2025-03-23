%% FIRI: Obstacle-Free Convex Polytope Generation Using Iterative MVIE
% This implementation computes restrictive halfspaces (RsI) and the Maximum Volume 
% Inscribed Ellipsoid (MVIE) for generating an obstacle-free convex region.
% The plot shows obstacles (light red), seed polytope (light green), seed point (blue), 
% halfspace boundaries (black dashed lines), and the MVIE ellipse (magenta).
%
% Make sure CasADi is installed and added to the MATLAB path.

import casadi.*;
clc; clear; close all;

%% Problem Setup
numObstacles = 5;
dim = 2; % 2D case
limits = [0 100; 0 100]; % [xmin xmax; ymin ymax]

% Generate non-overlapping obstacles (each obstacle is a convex polytope)
obstacles = generate_nonoverlapping_obstacles(numObstacles, dim, limits);

% Define seed: a point and its convex polytope Q (in V-representation)
seedPoint = [50, 50]; % Seed location
seedPolytope = [4, 4; 6, 4; 6, 6; 4, 6]; % Convex seed polytope Q

% Initialize ellipsoid E0 inside the seed region.
[AE, DE, bE] = initialize_ellipsoid(seedPolytope);

% Compute initial ellipsoid volume (area in 2D)
vol_prev = prod(diag(DE)) * (pi^(dim/2)/gamma(dim/2 + 1));
maxIter = 20;
tol = 1e-3;
converged = false;
iter = 0;

%% Main FIRI Iteration Loop
while ~converged && iter < maxIter
    iter = iter + 1;
    
    % RsI: Compute restrictive halfspaces for each obstacle in the transformed space.
    [AP, bP] = compute_restrictive_halfspaces(AE, DE, bE, obstacles, seedPolytope);
    
    % MVIE: Compute the Maximum Volume Inscribed Ellipsoid inside the polytope
    % defined by the selected halfspaces.
    [AE_new, DE_new, bE_new] = compute_mv_ie(AP, bP, dim, seedPolytope);
    
    % Compute current ellipsoid volume (area in 2D)
    vol_current = prod(diag(DE_new)) * (pi^(dim/2)/gamma(dim/2 + 1));
    rel_change = abs(vol_current - vol_prev) / vol_prev;
    
    % Check convergence: if relative change in volume is small enough, stop.
    if rel_change < tol
        converged = true;
    else
        AE = AE_new;
        DE = DE_new;
        bE = bE_new;
        vol_prev = vol_current;
    end
end

%% Final Halfspace Computation & Plotting of Results
[AP, bP] = compute_restrictive_halfspaces(AE, DE, bE, obstacles, seedPolytope);
plot_results(obstacles, seedPoint, seedPolytope, AP, bP, AE, DE, bE, limits);

%% Function Definitions

function obstacles = generate_nonoverlapping_obstacles(numObstacles, dim, limits)
    % Generate obstacles as convex polytopes that do not overlap.
    obstacles = cell(numObstacles, 1);
    min_dist = 2; 
    max_size = 1.5;
    centers = [];
    for i = 1:numObstacles
        while true
            new_center = rand(1, dim) .* (limits(:,2)' - limits(:,1)') + limits(:,1)';
            if isempty(centers) || min(vecnorm(centers - new_center, 2, 2)) > min_dist
                centers = [centers; new_center];
                break;
            end
        end
        % Generate points around the center and compute the convex hull.
        angles = linspace(0, 2*pi, 6)';
        radius = rand(6, 1) * max_size;
        pts = [radius .* cos(angles), radius .* sin(angles)] + new_center;
        k = convhull(pts);
        obstacles{i} = pts(k, :);
    end
end

function [AE, DE, bE] = initialize_ellipsoid(seedPolytope)
    % Initialize the ellipsoid E = { x = AE*DE*x0 + bE, ||x0||<=1 }
    center = mean(seedPolytope, 1);
    AE = eye(2);
    DE = diag([0.8, 0.8]);  % Initial semi-axis lengths (tunable)
    bE = center';
end

function [AP_selected, bP_selected] = compute_restrictive_halfspaces(AE, DE, bE, obstacles, seedPolytope)
    % Compute restrictive halfspaces by transforming the seed and obstacles into the
    % space where the current ellipsoid E becomes the unit ball.
    numObstacles = length(obstacles);
    dim = size(AE, 1);
    AP = zeros(numObstacles, dim);
    bP = zeros(numObstacles, 1);
    transformed_obstacles = cell(numObstacles, 1);
    
    % Transformation: x_transformed = (x - bE') / DE * AE'
    seed_transformed = (seedPolytope - bE') / DE * AE';
    for i = 1:numObstacles
        obs = obstacles{i};
        transformed_obstacles{i} = (obs - bE') / DE * AE';
    end
    
    % For each obstacle, solve:
    %   maximize ||a||^2  (i.e. minimize -||a||^2)
    % subject to:
    %   For each seed point:    a^T * x_seed <= a^T*a
    %   For each obstacle point: a^T * x_obs  >= a^T*a
    for i = 1:numObstacles
        opti = casadi.Opti();
        a = opti.variable(dim, 1);
        % Provide an initial guess for 'a' to help convergence.
        opti.set_initial(a, ones(dim,1));
        opti.minimize(-dot(a, a)); 
        
        % Seed inclusion constraints.
        for j = 1:size(seed_transformed, 1)
            opti.subject_to(seed_transformed(j, :) * a <= dot(a, a));
        end
        
        % Obstacle exclusion constraints.
        obs_transformed = transformed_obstacles{i};
        for j = 1:size(obs_transformed, 1)
            opti.subject_to(obs_transformed(j, :) * a >= dot(a, a));
        end
        
        % Set solver options.
        opti.solver('ipopt', struct(), struct('max_iter', 1000, 'tol', 1e-6, 'acceptable_tol', 1e-4, 'mu_init', 0.1, 'print_level', 0));
        try
            sol = opti.solve();
        catch ME
            warning('Solver did not converge for obstacle %d; using debug value.', i);
            sol = opti.debug;
        end
        a_opt = sol.value(a);
        AP(i, :) = a_opt';
        bP(i) = dot(a_opt, a_opt);  % bP = a^T*a defines the halfspace threshold
    end
    
    % Greedy selection of halfspaces.
    [~, order] = sort(sqrt(bP));
    active = true(numObstacles, 1);
    selected = [];
    
    for idx = 1:numObstacles
        i = order(idx);
        if ~active(i), continue; end
        
        selected(end+1) = i; %#ok<AGROW>
        a_i = AP(i, :)';
        thr = bP(i);
        
        for j = 1:numObstacles
            if active(j)
                outside = all(transformed_obstacles{j} * a_i >= thr - 1e-6);
                if outside
                    active(j) = false;
                end
            end
        end
        active(i) = false;
    end
    
    AP_selected = AP(selected, :);
    bP_selected = bP(selected);
end

function [AE, DE, bE] = compute_mv_ie(AP, bP, dim, seedPolytope)
    % Compute the Maximum Volume Inscribed Ellipsoid (MVIE) inside the convex polytope
    % defined by { x: AP*x <= bP }.
    opti = casadi.Opti();
    
    % Parameterize rotation in 2D via angle theta.
    theta = opti.variable();
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    
    % L: semi-axis lengths; c: ellipsoid center.
    L = opti.variable(dim, 1);
    c = opti.variable(dim, 1);
    
    % Objective: maximize the log-volume (in 2D, proportional to log(L1*L2)).
    opti.minimize(-sum(log(L + 1e-6)));
    
    % Enforce polytope containment for each halfspace.
    for i = 1:size(AP, 1)
        norm_term = norm(R * diag(L) * AP(i,:)');
        opti.subject_to(norm_term <= (bP(i) - AP(i,:)*c) - 1e-6);
    end
    
    % Ensure the ellipsoid is nondegenerate.
    opti.subject_to(L >= 1e-3);
    
    % Ensure the ellipsoid center is inside the polytope.
    opti.subject_to(AP*c <= bP - 1e-6);
    
    % Initialization using the centroid of the seed polytope.
    opti.set_initial(L, ones(dim,1)*0.5);
    opti.set_initial(c, mean(seedPolytope)');  
    opti.set_initial(theta, 0);
    
    % Set solver options.
    opts = struct('max_iter', 1000, 'tol', 1e-6, 'acceptable_tol', 1e-4, 'mu_init', 0.1, 'print_level', 0);
    opti.solver('ipopt', struct(), opts);
    
    sol = opti.solve();
    
    % Extract solution.
    theta_opt = sol.value(theta);
    AE = [cos(theta_opt) -sin(theta_opt); sin(theta_opt) cos(theta_opt)];
    DE = diag(sol.value(L));
    bE = sol.value(c);
end

function plot_results(obstacles, seedPoint, seedPolytope, AP, bP, AE, DE, bE, limits)
    % Plot obstacles, seed, halfspaces, and the MVIE.
    figure; hold on;
    xlim(limits(1,:)); ylim(limits(2,:));
    axis equal; grid on;
    title('FIRI: Obstacle-Free Convex Space Generation');
    xlabel('X'); ylabel('Y');
    
    % Plot obstacles in light red.
    for i = 1:length(obstacles)
        fill(obstacles{i}(:,1), obstacles{i}(:,2), [1 0.6 0.6], 'EdgeColor', 'r', 'FaceAlpha', 0.5);
    end
    
    % Plot seed polytope in light green.
    fill(seedPolytope(:,1), seedPolytope(:,2), [0.6 1 0.6], 'EdgeColor', 'g', 'FaceAlpha', 0.5);
    
    % Plot seed point in blue.
    scatter(seedPoint(1), seedPoint(2), 100, 'b', 'filled');
    
    % Plot each halfspace boundary.
    x_vals = linspace(limits(1,1), limits(1,2), 100);
    for i = 1:size(AP,1)
        % Check if the second component is nearly zero.
        if abs(AP(i,2)) < 1e-3
            % For a nearly vertical halfspace, plot a vertical line.
            x_line = repmat(bP(i)/AP(i,1), size(x_vals));
            plot(x_line, linspace(limits(2,1), limits(2,2), 100), 'k--', 'LineWidth', 1.5);
        else
            % Otherwise, use the line equation: AP(i,1)*x + AP(i,2)*y = bP(i)
            y_vals = (bP(i) - AP(i,1)*x_vals) / AP(i,2);
            plot(x_vals, y_vals, 'k--', 'LineWidth', 1.5);
        end
    end
    
    % Plot the final MVIE ellipse in magenta.
    theta_plot = linspace(0, 2*pi, 200);
    ellipse_pts = AE * DE * [cos(theta_plot); sin(theta_plot)] + bE;
    plot(ellipse_pts(1,:), ellipse_pts(2,:), 'm-', 'LineWidth', 2);
    
    legend({'Obstacles', 'Seed Polytope', 'Seed Point', 'Halfspaces', 'MVIE'}, 'Location', 'best');
    hold off;
end
