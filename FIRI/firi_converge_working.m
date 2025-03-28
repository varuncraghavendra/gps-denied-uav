function FIRI_main
    % FIRI_main.m
    % This script implements a 2-D version of the FIRI algorithm (Algorithm 1)
    % for computing an obstacle-free convex polytope and its maximum-volume
    % inscribed ellipse. It uses CasADi for the MVIE optimization and includes
    % plotting of:
    %   (i) non-overlapping convex obstacles,
    %   (ii) a seed ellipsoid (and its convex hull),
    %   (iii) the final inflated polytope, and
    %   (iv) the final inscribed ellipse.
    %
    % The MVIE is computed by parameterizing the ellipse as
    %   E = { x : norm(L^{-1}(x-d)) <= 1 }
    % with L lower-triangular (L(1,2)==0) and d in R^2.
    % The code enforces that the halfspace constraints computed in the RsI
    % module (via tangent-plane computation) do not intrude into the obstacles.
    
    close all; clear; clc;
    import casadi.*
    
    %% Environment parameters (in meters)
    env_size = [50, 50];          % 50m x 50m environment
    num_obstacles = 30;           % Number of convex obstacles
    min_spacing = 0.3;            % Minimum spacing between obstacles
    ellipsoid_spacing = 2.0;      % Minimum spacing from the ellipsoid seed
    
    %% Ellipsoid (seed) parameters
    % For 2-D, we set AE = 1 (no rotation) and specify the semi-axis lengths.
    AE = 1;                     
    DE = diag([4, 2]);          % Semi-axis lengths: 4m and 2m
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
    % Initial ellipse for FIRI: a very small ellipse around the seed,
    % represented by L = 0.01*eye(2) and center = seed.
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
    draw_ellipse(E_final.L, E_final.d);  % Plot final inscribed ellipse.
    
    hold off;
    
    
    %% ---------------- Nested Functions ---------------- %%
    
    %% FIRI: Outer loop alternating between RsI and MVIE.
    function [Pk, E_final, iter] = FIRI(seed_point, obstacles, L0, d0, tol, max_iters)
        % The ellipse is parameterized as: E = { x : norm(L^{-1}(x-d)) <= 1 }.
        % Initialize with L0 and d0.
        E.L = L0;  E.d = d0;
        iter = 0;
        % For 2-D, use the product of the diagonal entries as a proxy for area.
        vol_prev = E.L(1,1) * E.L(2,2);
        
        while iter < max_iters
            iter = iter + 1;
            % ----- RsI Module -----
            [A, b] = RsI(E, obstacles, seed_point);
            % Define polytope Pk = { x : A*x <= b }.
            Pk = struct('A', A, 'b', b);
            
            % Check that the seed is contained.
            if ~isempty(A) && any(A' * seed_point >= b)
                disp('Terminating early: seed point not contained in region.');
                break;
            end
            
            % ----- MVIE Module -----
            [L_new, d_new] = MVIE(Pk);
            E_new.L = L_new;  E_new.d = d_new;
            vol_new = L_new(1,1) * L_new(2,2);
            disp(['Iteration ', num2str(iter), ', Volume proxy = ', num2str(vol_new)]);
            
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
    % For each obstacle, the function computes the boundary point that is
    % closest to the seed (using edge-projection) and then constructs a tangent plane.
    function [A, b] = RsI(E, obstacles, seed_point)
        A_list = [];
        b_list = [];
        for i = 1:length(obstacles)
            obs = obstacles{i};
            if size(unique(obs, 'rows'), 1) < 3, continue; end
            % Compute the closest point on the obstacle boundary to the seed.
            x_star = ClosestPointOnPolygon(seed_point, obs);
            % Compute the tangent plane at x_star.
            [a, bi] = TangentPlane(seed_point, x_star);
            A_list = [A_list, a];
            b_list = [b_list; bi];
        end
        A = A_list;
        b = b_list;
    end

    %% ClosestPointOnPolygon: For a convex polygon, compute the boundary point
    % that minimizes the Euclidean distance to the seed.
    function x_star = ClosestPointOnPolygon(seed, poly)
        num_pts = size(poly,1);
        min_dist = inf;
        x_star = poly(1,:)';
        for j = 1:num_pts
            p1 = poly(j,:)';
            p2 = poly(mod(j, num_pts)+1,:)';
            % Project seed onto the edge (p1,p2)
            v = p2 - p1;
            if norm(v) < eps, continue; end
            t = max(0, min(1, (seed - p1)' * v / (v' * v)));
            proj = p1 + t*v;
            d = norm(seed - proj);
            if d < min_dist
                min_dist = d;
                x_star = proj;
            end
        end
    end

    %% TangentPlane: Compute the tangent plane (halfspace) at the contact point.
    % The halfspace is defined as { x : a'*x <= a'*x_star },
    % where the normal vector is computed as (x_star - seed).
    function [a, bi] = TangentPlane(seed, x_star)
        a = x_star - seed;
        % Scale a such that the halfspace is tight.
        % Here, we do not normalize a so that a'*x_star is the correct offset.
        bi = a' * x_star;
    end

    %% MVIE: Compute the maximum-volume inscribed ellipse via an SOCP formulation using CasADi.
    % The ellipse is parameterized as: E = { x : norm(L^{-1}(x-d)) <= 1 },
    % with L lower-triangular (L(1,2)==0).
    function [L_new, d_new] = MVIE(P)
        import casadi.*
        n = 2;
        % Decision variables: L (2x2) and d (2x1).
        L_sym = SX.sym('L', n, n);
        d_sym = SX.sym('d', n, 1);
        % Impose lower-triangular structure: L(1,2)==0.
        constr_L = L_sym(1,2);
        % Objective: maximize area ~ log(L(1,1))+log(L(2,2)) (minimize negative).
        obj = - ( log(L_sym(1,1)) + log(L_sym(2,2)) );
        
        g = [];
        A = P.A;
        b_vec = P.b;
        num_con = length(b_vec);
        for i = 1:num_con
            a_i = A(:,i);
            % Constraint: norm(L_sym'*a_i,2) + a_i' * d_sym <= b_vec(i)
            g = [g; norm(L_sym'*a_i,2) + a_i' * d_sym - b_vec(i)];
        end
        % Enforce lower-triangular structure: L(1,2)==0.
        g = [g; constr_L];
        
        x_dec = [reshape(L_sym, n*n, 1); d_sym];
        nlp.x = x_dec;
        nlp.f = obj;
        nlp.g = g;
        opts.ipopt.print_level = 0; opts.print_time = false;
        solver = nlpsol('solver','ipopt', nlp, opts);
        % Initial guess: L = 0.1*eye(2), d = seed.
        x0 = [reshape(eye(n)*0.1, n*n, 1); seed_point];
        sol = solver('x0', x0, 'lbg', -inf*ones(length(g),1), 'ubg', zeros(length(g),1));
        solx = full(sol.x);
        L_new = reshape(solx(1:n*n), n, n);
        d_new = solx(n*n+1:end);
    end

    %% hrep2vrep: Convert H-representation (A*x <= b) to vertices (2-D).
    function vertices = hrep2vrep(A, b)
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
            center = mean(pts,1);
            angles = atan2(pts(:,2)-center(2), pts(:,1)-center(1));
            [~, order] = sort(angles);
            vertices = pts(order,:);
        else
            vertices = [];
        end
    end

    %% draw_ellipse: Plot ellipse defined as E = { x : norm(L^{-1}(x-d)) <= 1 }.
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
        Q = (AE * DE * x)';                % 100 x 2
        Q = Q + bE';                       % translate by bE
    end

    %% generate_nonoverlapping_obstacles: Generate convex obstacles using convex hulls.
    function convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_obs, min_spacing, Q, ellipsoid_spacing)
        convex_obstacles = cell(num_obs, 1);
        placed_centers = zeros(num_obs, 2);
        for i = 1:num_obs
            num_pts = randi([3, 8]);
            is_valid = false;
            while ~is_valid
                center = [rand * env_size(1), rand * env_size(2)];
                if (i == 1 || all(vecnorm(placed_centers(1:i-1,:) - center, 2, 2) > min_spacing)) && ...
                   all(vecnorm(Q - center, 2, 2) > ellipsoid_spacing)
                    is_valid = true;
                    placed_centers(i, :) = center;
                    angles = linspace(0, 2*pi, num_pts+1)' + rand*pi/4;
                    angles(end) = [];
                    radii = rand(num_pts, 1) * 2 + 1;
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
