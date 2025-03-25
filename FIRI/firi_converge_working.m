function FIRI_main
    % FIRI_main.m
    % This MATLAB script implements a 2-D version of FIRI (Algorithm 1)
    % for generating an obstacle-free convex polytope and computing its
    % maximum-volume inscribed ellipse using CasADi for optimization.
    %
    % Save this file as FIRI_main.m and run FIRI_main in MATLAB.
    % Requirements:
    % - CasADi must be installed and added to your MATLAB path.
    % - A function con2vert is assumed to exist for converting an H-
    %   representation to vertices (or modify hrep2vrep accordingly).

    close all; clear; clc;
    import casadi.*
    
    %% Environment and initial parameters
    obstacles = gen_obstacles();
    seed_point = [0.5; 0.5];  % seed point in [0,1]x[0,1]
    % Initial ellipse: very small ellipse around the seed.
    C0 = eye(2)*0.01;  
    d0 = seed_point;
    tol = 1e-5;
    max_iters = 10;

    %% Plot initial obstacles and seed
    figure; hold on; grid on; axis equal;
    xlim([0 1]); ylim([0 1]);
    title('FIRI: Obstacle Environment and Final Polytope');
    xlabel('X'); ylabel('Y');
    plot_obstacles(obstacles);
    plot(seed_point(1), seed_point(2), 'ko', 'MarkerFaceColor','k');
    
    %% Run FIRI algorithm
    [Pk, E_final, iter] = FIRI(seed_point, obstacles, C0, d0, tol, max_iters);
    disp(['FIRI converged in ' num2str(iter) ' iterations.']);
    
    %% Plot final polytope and ellipse
    if ~isempty(Pk)
        plot_polytope(Pk);
    end
    draw_ellipse(E_final.C, E_final.d);
    
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIRI: Alternating between RsI and MVIE
    function [Pk, E_final, iter] = FIRI(seed_point, obstacles, C0, d0, tol, max_iters)
        import casadi.*
        % E is a struct representing the ellipse:
        % E = { x | norm(C^{-1}(x-d)) <= 1 }
        E.C = C0;
        E.d = d0;
        iter = 0;
        vol_prev = det(E.C);  % Using det(C) as a proxy for area
        
        while iter < max_iters
            iter = iter + 1;
            % RsI: Compute halfspace constraints from obstacles using the current ellipse E.
            [A, b] = RsI(E, obstacles, seed_point);
            % Polytope P is given in H-representation: { x | A'*x <= b }.
            Pk = struct('A', A, 'b', b);
            
            % Check seed containment only if A is non-empty.
            if ~isempty(A) && any(A' * seed_point >= b)
                disp('Terminating early: seed point not contained in region.');
                break;
            end
            
            % MVIE: Compute maximum-volume inscribed ellipse in polytope Pk using CasADi.
            [C_new, d_new] = MVIE(Pk);
            E_new.C = C_new; E_new.d = d_new;
            vol_new = det(E_new.C);
            disp(['Iteration ', num2str(iter) ', det(C) = ', num2str(vol_new)]);
            
            if (vol_new - vol_prev)/vol_prev < tol
                E_final = E_new;
                break;
            end
            
            E = E_new;
            vol_prev = vol_new;
        end
        E_final = E;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RsI: Compute separating halfspaces from obstacles.
    function [A, b] = RsI(E, obstacles, seed_point)
        A_list = [];
        b_list = [];
        % Loop over each obstacle (each represented as a set of vertices)
        for i = 1:length(obstacles)
            obs = obstacles{i};  % obs is an N_i x 2 array
            if size(unique(obs, 'rows'), 1) < 3
                % Skip degenerate obstacles.
                continue;
            end
            [x_star, ~] = ClosestPointOnObstacle(E, obs);
            [a, bi] = TangentPlane(E, x_star);
            A_list = [A_list, a];
            b_list = [b_list; bi];
        end
        A = A_list;
        b = b_list;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ClosestPointOnObstacle: Solve a QP via CasADi to compute the closest point.
    function [x_star, dist] = ClosestPointOnObstacle(E, obs)
        % Transform obstacle vertices into the ellipse's coordinate space.
        C_inv = inv(E.C);
        % Subtract E.d from each column of obs'
        obs_centered = obs' - repmat(E.d, 1, size(obs,1));
        obs_tilde = C_inv * obs_centered;  % 2 x m matrix
        obs_tilde = obs_tilde';  % m x 2: each row is a transformed vertex.
        [m, ~] = size(obs_tilde);
        
        import casadi.*
        w = SX.sym('w', m, 1);  % Weight variables for convex combination.
        x_tilde = obs_tilde' * w;  % 2 x 1 vector.
        objective = dot(x_tilde, x_tilde);
        
        % Constraints: sum(w) == 1 and w >= 0.
        g = [sum1(w) - 1; w];
        
        % Set up the NLP (a QP).
        nlp.x = w;
        nlp.f = objective;
        nlp.g = g;
        opts.ipopt.print_level = 0;
        opts.print_time = false;
        solver = nlpsol('solver', 'ipopt', nlp, opts);
        % Initial guess:
        w0 = ones(m,1)/m;
        sol = solver('x0', w0, 'lbg', [1; zeros(m,1)], 'ubg', [1; Inf(m,1)]);
        w_opt = full(sol.x);
        obj_val = full(sol.f);
        
        % Distance is sqrt(objective) - 1 (since unit ball boundary is 1).
        dist = sqrt(obj_val) - 1;
        x_star = E.C * (obs_tilde' * w_opt) + E.d;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TangentPlane: Compute halfspace from the contact point.
    function [a, bi] = TangentPlane(E, x_star)
        C_inv = inv(E.C);
        C_inv2 = C_inv * C_inv';
        a = 2 * C_inv2 * (x_star - E.d);
        bi = a' * x_star;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MVIE: Compute maximum-volume inscribed ellipse via SOCP formulation.
    function [C_new, d_new] = MVIE(P)
        % We solve: maximize log(det(C)) subject to:
        %   norm(a_i'*C,2) + a_i'*d <= b_i, for each halfspace constraint.
        % Decision variables: C (2x2 matrix) and d (2x1 vector).
        import casadi.*
        n = 2;
        C_var = SX.sym('C', n, n);
        d_var = SX.sym('d', n, 1);
        
        % Objective: maximize log(det(C)) (we minimize -log(det(C))).
        obj = -log(det(C_var));
        
        g = [];
        A = P.A;
        b = P.b;
        num_con = length(b);
        for i = 1:num_con
            a_i = A(:,i);
            % Constraint: norm(a_i'*C_var,2) + a_i'*d_var <= b(i).
            g = [g; norm(a_i' * C_var, 2) + a_i' * d_var - b(i)];
        end
        
        % Build the NLP.
        x = [reshape(C_var, n*n, 1); d_var];
        nlp.x = x;
        nlp.f = obj;
        nlp.g = g;
        opts.ipopt.print_level = 0;
        opts.print_time = false;
        solver = nlpsol('solver', 'ipopt', nlp, opts);
        
        % Initial guess: a small ellipse near the origin.
        x0 = [reshape(eye(n)*0.1, n*n, 1); zeros(n,1)];
        sol = solver('x0', x0, 'lbg', -Inf*ones(length(g),1), 'ubg', zeros(length(g),1));
        solx = full(sol.x);
        C_new = reshape(solx(1:n*n), n, n);
        d_new = solx(n*n+1:end);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Obstacle generation using random points and MATLAB's alphashape.
    function obstacles = gen_obstacles()
        n_points = 80;
        alpha = 15;
        points = rand(n_points, 2);
        % Create an alphaShape. Adjust the parameter as needed.
        shp = alphaShape(points, 1/alpha);
        [tri, pts] = boundaryFacets(shp);
        obstacles = {};
        for i = 1:size(tri,1)
            obs_i = pts(tri(i,:), :);
            % Only add obstacles with at least 3 unique points.
            if size(unique(obs_i, 'rows'), 1) >= 3
                obstacles{end+1} = obs_i;
            end
        end
        % Replicate obstacles shifted in various directions to mimic periodicity.
        orig = obstacles;
        shifts = [1,0; -1,0; 0,1; 0,-1; 1,1; 1,-1; -1,1; -1,-1];
        for s = 1:size(shifts,1)
            for i = 1:length(orig)
                obs_shifted = orig{i} + shifts(s,:);
                if size(unique(obs_shifted, 'rows'), 1) >= 3
                    obstacles{end+1} = obs_shifted;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting functions
    function plot_obstacles(obstacles)
        for i = 1:length(obstacles)
            obs = obstacles{i};
            % Only plot if there are at least three unique points.
            if size(unique(obs, 'rows'), 1) >= 3
                try
                    k = convhull(obs(:,1), obs(:,2));
                    fill(obs(k,1), obs(k,2), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
                catch ME
                    warning('Skipping obstacle due to convhull error: %s', ME.message);
                end
            end
        end
    end

    function plot_polytope(P)
        % Convert H-representation (A'*x <= b) into vertices.
        vertices = hrep2vrep(P.A', P.b);
        if ~isempty(vertices)
            try
                k = convhull(vertices(:,1), vertices(:,2));
                plot(vertices(k,1), vertices(k,2), 'b-', 'LineWidth', 2);
            catch ME
                warning('Error computing convex hull for polytope: %s', ME.message);
            end
        else
            warning('Polytope vertices could not be computed.');
        end
    end

    function draw_ellipse(C, d)
        theta = linspace(0, 2*pi, 100);
        circle = [cos(theta); sin(theta)];
        ellipse = C * circle + d;
        plot(ellipse(1,:), ellipse(2,:), 'g-', 'LineWidth', 2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% hrep2vrep: Convert H-representation (A'*x <= b) to vertices.
    % This implementation uses the external function con2vert if available.
    function vertices = hrep2vrep(A, b)
        % A: m x n, b: m x 1; polytope { x in R^n: A*x <= b }
        try
            vertices = con2vert([A, b]);
        catch
            vertices = [];
            warning('con2vert function not found. Please supply your own hrep2vrep implementation.');
        end
    end

end
