function FIRI_main
    % FIRI_main.m
    % This script implements a 2-D version of the FIRI algorithm (Algorithm 1)
    % for computing an obstacle-free convex polytope and its maximum-volume
    % inscribed ellipse. It uses CasADi for optimization and includes plotting
    % of:
    %   (i) non-overlapping convex obstacles,
    %   (ii) a seed ellipsoid (and its convex hull),
    %   (iii) the final inflated polytope, and
    %   (iv) the final inscribed ellipse.
    %
    % The MVIE is computed by parameterizing the ellipse as
    %   E = { x : norm(L^{-1}(x-d)) <= 1 }
    % with L lower-triangular (L(1,2)==0) and d in R^2. The objective is to
    % maximize log(L(1,1))+log(L(2,2)) subject to the constraints
    %   norm(a_i'*L,2) + a_i'*d <= b_i, for each halfspace from the obstacles.
    %
    % Save this file as FIRI_main.m and run it in MATLAB.
    
    close all; clear; clc;
    import casadi.*
    
    %% Environment parameters (in meters)
    env_size = [50, 50];          % 50m x 50m environment
    num_obstacles = 30;           % Number of convex obstacles
    min_spacing = 0.3;            % Minimum spacing between obstacles
    ellipsoid_spacing = 2.0;      % Minimum spacing from the ellipsoid to obstacles
    
    %% Ellipsoid (seed) parameters
    AE = 1;                     % For 2-D, use scalar 1 for rotation
    DE = diag([4, 2]);          % Semi-axis lengths (4m and 2m)
    bE = [15; 25];              % Center of the seed ellipsoid at (15,25)
    
    % Generate seed ellipsoid points Q and its convex hull P_seed
    Q = generate_ellipsoid_seed(AE, DE, bE);
    P_seed = convhull(Q(:,1), Q(:,2));
    
    %% Generate non-overlapping convex obstacles avoiding the ellipsoid seed Q
    convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_obstacles, min_spacing, Q, ellipsoid_spacing);
    
    %% Plot the environment: obstacles, seed ellipsoid and its convex hull
    figure; hold on; grid on; axis equal;
    xlim([0, env_size(1)]); ylim([0, env_size(2)]);
    xlabel('X (m)'); ylabel('Y (m)');
    title('Environment: Convex Obstacles, Seed Ellipsoid & Inflated Region');
    
    % Plot obstacles (red)
    for i = 1:numel(convex_obstacles)
        obs = convex_obstacles{i};
        if size(unique(obs, 'rows'), 1) >= 3
            k = convhull(obs(:,1), obs(:,2));
            fill(obs(k,1), obs(k,2), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
        end
    end
    
    % Plot convex hull of seed ellipsoid (blue)
    plot(Q(P_seed,1), Q(P_seed,2), 'b-', 'LineWidth', 2);
    % Plot seed ellipsoid points (green circles)
    plot(Q(:,1), Q(:,2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    
    %% Run FIRI algorithm using the seed point bE.
    seed_point = bE;
    % For FIRI, use the convex obstacles generated above.
    % Initial ellipse for FIRI: very small ellipse around seed (represented by L = 0.1*eye(2) and center = seed)
    [Pk, E_final, iter] = FIRI(seed_point, convex_obstacles, eye(2)*0.01, seed_point, 1e-5, 10);
    disp(['FIRI converged in ' num2str(iter) ' iterations.']);
    
    %% Plot final inflated polytope and final inscribed ellipse.
    if ~isempty(Pk)
        vertices = hrep2vrep(Pk.A', Pk.b);
        if ~isempty(vertices)
            k = convhull(vertices(:,1), vertices(:,2));
            plot(vertices(k,1), vertices(k,2), 'm-', 'LineWidth', 2);
        else
            warning('Final polytope vertices could not be computed.');
        end
    end
    draw_ellipse(E_final.L, E_final.d);  % Here E_final.L is the L factor
    
    hold off;
    
    %% ---------------- Nested Functions ---------------- %%
    
    %% FIRI: Alternates between Restrictive Inflation (RsI) and MVIE.
    function [Pk, E_final, iter] = FIRI(seed_point, obstacles, L0, d0, tol, max_iters)
        % Here the ellipse is parameterized as E = { x : norm(L^{-1}(x-d)) <= 1 }.
        % We initialize with L0 and d0.
        E.L = L0;
        E.d = d0;
        iter = 0;
        vol_prev = E.L(1,1) * E.L(2,2);  % Use product of diagonal entries as a proxy for area
        
        while iter < max_iters
            iter = iter + 1;
            % RsI: Compute halfspace constraints from obstacles given current ellipse E.
            [A, b] = RsI(E, obstacles, seed_point);
            % Polytope P in H-representation: { x : A*x <= b }.
            % Note: In our code, we use the constraints in the form
            %       a_i'*x <= b_i, where a_i comes from the tangent plane.
            Pk = struct('A', A, 'b', b);
            
            % Check seed containment if any constraint is violated.
            if ~isempty(A) && any(A' * seed_point >= b)
                disp('Terminating early: seed point not contained in region.');
                break;
            end
            
            % MVIE: Compute maximum-volume inscribed ellipse in polytope Pk.
            [L_new, d_new] = MVIE(Pk);
            E_new.L = L_new; E_new.d = d_new;
            vol_new = L_new(1,1) * L_new(2,2);
            disp(['Iteration ', num2str(iter) ', Volume proxy = ', num2str(vol_new)]);
            
            if (vol_new - vol_prev) / vol_prev < tol
                E_final = E_new;
                break;
            end
            
            E = E_new;
            vol_prev = vol_new;
        end
        E_final = E;
    end

    %% RsI: Compute separating halfspaces from obstacles.
    function [A, b] = RsI(E, obstacles, seed_point)
        A_list = [];
        b_list = [];
        % Loop over each obstacle (each is a set of vertices)
        for i = 1:length(obstacles)
            obs = obstacles{i};
            if size(unique(obs, 'rows'), 1) < 3
                continue;  % Skip degenerate obstacles.
            end
            [x_star, ~] = ClosestPointOnObstacle(E, obs);
            [a, bi] = TangentPlane(E, x_star);
            A_list = [A_list, a];
            b_list = [b_list; bi];
        end
        A = A_list;
        b = b_list;
    end

    %% ClosestPointOnObstacle: Solve a QP via CasADi to compute the closest point on an obstacle.
    function [x_star, dist] = ClosestPointOnObstacle(E, obs)
        % Transform obstacle vertices into the ellipse coordinate space.
        C = E.L; % Here we use L as our scale factor (since our ellipse is defined by L)
        % To transform, note that if E = { x : norm(L^{-1}(x-d)) <= 1 },
        % then x -> L^{-1}(x-d). Compute:
        L_inv = inv(E.L);
        obs_centered = obs' - repmat(E.d, 1, size(obs,1));
        obs_tilde = L_inv * obs_centered;  % 2 x m
        obs_tilde = obs_tilde';             % m x 2
        [m, ~] = size(obs_tilde);
        
        import casadi.*
        w = SX.sym('w', m, 1);
        x_tilde = obs_tilde' * w;  % 2 x 1
        objective = dot(x_tilde, x_tilde);
        
        % Constraints: sum(w)==1, w >= 0.
        g = [sum1(w) - 1; w];
        nlp.x = w;
        nlp.f = objective;
        nlp.g = g;
        opts.ipopt.print_level = 0;
        opts.print_time = false;
        solver = nlpsol('solver', 'ipopt', nlp, opts);
        w0 = ones(m,1)/m;
        sol = solver('x0', w0, 'lbg', [1; zeros(m,1)], 'ubg', [1; Inf(m,1)]);
        w_opt = full(sol.x);
        obj_val = full(sol.f);
        dist = sqrt(obj_val) - 1;
        x_star = E.L * (obs_tilde' * w_opt) + E.d;
    end

    %% TangentPlane: Compute halfspace from the contact point.
    function [a, bi] = TangentPlane(E, x_star)
        % For our ellipse parameterization, we compute the tangent plane as:
        % a = 2 * (x_star - d) scaled by an appropriate inverse metric.
        % Here we use a simple formulation similar to the CVX version.
        L_inv = inv(E.L);
        L_inv2 = L_inv * L_inv';
        a = 2 * L_inv2 * (x_star - E.d);
        bi = a' * x_star;
    end

    %% MVIE: Compute the maximum-volume inscribed ellipse via an SOCP formulation.
    % This function is reparameterized so that the ellipse is given by:
    %    E = { x : norm(L^{-1}(x-d)) <= 1 }
    % where L is a 2x2 lower-triangular matrix (with L(1,2)==0).
    function [L_new, d_new] = MVIE(P)
        import casadi.*
        n = 2;
        % Decision variables: L (2x2) and d (2x1).
        L_var = SX.sym('L', n, n);
        d_var = SX.sym('d', n, 1);
        % Impose lower-triangular structure: L(1,2)==0.
        constr_L = L_var(1,2);
        % Objective: maximize volume ~ log(L(1,1))+log(L(2,2)) i.e. minimize negative.
        obj = - (log(L_var(1,1)) + log(L_var(2,2)));
        
        g = [];
        A = P.A;
        b = P.b;
        num_con = length(b);
        for i = 1:num_con
            a_i = A(:,i);
            r = a_i' * L_var;  % 1x2 row vector
            % Constraint: norm(r,2) + a_i' * d_var <= b(i)
            g = [g; norm(r,2) + a_i' * d_var - b(i)];
        end
        % Add constraint that L(1,2)==0.
        g = [g; constr_L];
        
        x = [reshape(L_var, n*n, 1); d_var];
        nlp.x = x;
        nlp.f = obj;
        nlp.g = g;
        opts.ipopt.print_level = 0;
        opts.print_time = false;
        solver = nlpsol('solver','ipopt', nlp, opts);
        % Initial guess: small ellipse: L = 0.1*eye(2), d = seed_point.
        x0 = [reshape(eye(n)*0.1, n*n, 1); seed_point];
        sol = solver('x0', x0, 'lbg', -Inf*ones(length(g),1), 'ubg', zeros(length(g),1));
        solx = full(sol.x);
        L_new = reshape(solx(1:n*n), n, n);
        d_new = solx(n*n+1:end);
    end

    %% hrep2vrep: Convert H-representation (A*x <= b) to vertices (2-D version).
    function vertices = hrep2vrep(A, b)
        % Here A is m x 2 and b is m x 1; we compute all intersections.
        tol = 1e-6;
        m = size(A,1);
        pts = [];
        for i = 1:m-1
            for j = i+1:m
                M = [A(i,:); A(j,:)];
                if rank(M) < 2, continue; end
                x_int = M \ [b(i); b(j)];
                if all(A * x_int <= b + tol)
                    pts = [pts; x_int'];
                end
            end
        end
        if ~isempty(pts)
            pts = unique(pts, 'rows');
            centroid = mean(pts,1);
            angles = atan2(pts(:,2)-centroid(2), pts(:,1)-centroid(1));
            [~, order] = sort(angles);
            vertices = pts(order,:);
        else
            vertices = [];
        end
    end

    %% draw_ellipse: Plot ellipse defined by L and d.
    function draw_ellipse(L, d)
        theta = linspace(0, 2*pi, 100);
        circle = [cos(theta); sin(theta)];
        ellipse = L * circle + d;
        plot(ellipse(1,:), ellipse(2,:), 'g-', 'LineWidth', 2);
    end

    %% generate_ellipsoid_seed: Generate 100 points on the ellipse defined by AE, DE, and bE.
    function Q = generate_ellipsoid_seed(AE, DE, bE)
        theta = linspace(0, 2*pi, 100)';  % 100 points on unit circle
        x = [cos(theta), sin(theta)]';     % 2 x 100
        Q = (AE * DE * x)';                % 100 x 2 (scale)
        Q = Q + bE';                       % translate by bE
    end

    %% generate_nonoverlapping_obstacles: Generate convex obstacles using convex hulls.
    function convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_obstacles, min_spacing, Q, ellipsoid_spacing)
        convex_obstacles = cell(num_obstacles, 1);
        placed_centers = zeros(num_obstacles, 2);
        for i = 1:num_obstacles
            num_points = randi([3, 8]);
            is_valid = false;
            while ~is_valid
                center = [rand * env_size(1), rand * env_size(2)];
                if (i == 1 || all(vecnorm(placed_centers(1:i-1,:) - center, 2, 2) > min_spacing)) && ...
                   all(vecnorm(Q - center, 2, 2) > ellipsoid_spacing)
                    is_valid = true;
                    placed_centers(i, :) = center;
                    angles = linspace(0, 2*pi, num_points+1)' + rand*pi/4;
                    angles(end) = [];
                    radii = rand(num_points, 1) * 2 + 1;
                    x_obs = center(1) + radii .* cos(angles);
                    y_obs = center(2) + radii .* sin(angles);
                    x_obs = min(max(x_obs, 0), env_size(1));
                    y_obs = min(max(y_obs, 0), env_size(2));
                    convex_obstacles{i} = [x_obs, y_obs];
                end
            end
        end
    end

end
