function Pk = firi_ell(A, b, lb, ub, Q, O, E0, rho)
    % firi_ell: Fast Iterative Region Inflation for ellipsoidal regions.
    %
    % Inputs:
    %   A, b   - Constraints for the region of interest.
    %   lb, ub - Lower and upper bounds for the region.
    %   Q      - Seed points. Expected as an n x m matrix (n = space dimension)
    %            or as an m x n matrix.
    %   O      - Obstacles (cell array; each cell is an mₒ x n matrix, one obstacle).
    %   E0     - (Optional) Initial ellipsoid structure.
    %   rho    - (Optional) Convergence threshold (default: 0.01).
    %
    % Output:
    %   Pk     - Convex polytope computed by the algorithm.
    
    % Determine space dimensionality from Q.
    % Accept Q as either n x m or m x n.
    if size(Q,2) < size(Q,1)
        Q = Q';  % Transpose if needed.
    end
    n = size(Q,1);  % Now Q is n x m.
    
    % Default values for E0 and rho if not provided.
    if nargin < 7 || isempty(E0)
        E0 = initialize_ellipsoid(n);
    end
    if nargin < 8 || isempty(rho)
        rho = 0.01;
    end

    % Save the initial ellipsoid for plotting.
    E_initial = E0;
    
    % Initialize variables.
    k = 0;
    E = E0;         % Starting ellipsoid.
    P = [];         % Placeholder for convex polytope (in halfspace form).

    % Main iterative loop.
    while true
        k = k + 1;
        
        % RsI step: transform the seed and obstacles into the ellipsoid's space.
        Q_transformed = transform_seed(E, Q);      % Q is now n x m.
        O_transformed = transform_obstacles(E, O);   % Each obstacle: mₒ x n.
        
        % For each obstacle, compute a halfspace.
        for i = 1:length(O_transformed)
            % Use CasADi to maximize the halfspace direction.
            ai = maximize_halfspace(Q_transformed, O_transformed{i});
            % Add the halfspace (represented by ai) to the polytope.
            P = add_halfspace(P, ai);
        end
        
        % Update polytope Pk (recover to original space).
        Pk = update_polytope(P, E);
        
        % Compute the Maximum Volume Inscribed Ellipsoid (MVIE) of Pk.
        Ek = compute_mvie(Pk);
        
        % Check convergence based on ellipsoid volume.
        if volume(Ek) <= (1 + rho) * volume(E)
            break;
        end
        
        E = Ek;  % Update ellipsoid for next iteration.
    end
    
    % Plot the results (for 2D only)
    if n == 2
        plot_firi_results(E_initial, Q, O, Pk);
    else
        disp('Plotting is implemented only for 2D.');
    end
end

%% Helper Functions

% Initialize a random ellipsoid in n-dimensional space.
function E = initialize_ellipsoid(n)
    AE = eye(n);                   % Identity as the orthonormal matrix.
    DE = diag(rand(n, 1) + 0.5);     % Diagonal with positive entries.
    bE = rand(n, 1) * 10;            % Random translation vector.
    E = struct('AE', AE, 'DE', DE, 'bE', bE);
end

% Transform seed Q using ellipsoid parameters.
function Q_transformed = transform_seed(E, Q)
    % Q is expected to be n x m (n dimensions, m points).
    n = size(E.AE, 1);
    if size(Q,1) ~= n
        if size(Q,2) == n
            Q = Q';
        else
            error('Seed dimensions do not match ellipsoid dimensionality.');
        end
    end
    % Subtract bE from every column.
    Q_transformed = (inv(E.DE) * E.AE') * (Q - repmat(E.bE, 1, size(Q,2)));
end

% Transform each obstacle in O using the ellipsoid parameters.
function O_transformed = transform_obstacles(E, O)
    % Each obstacle is expected to be an mₒ x n matrix (each row is an n-dimensional point).
    n = size(E.AE, 1);
    O_transformed = cell(size(O));
    for i = 1:length(O)
        % If an obstacle has its points in columns, transpose it.
        if size(O{i}, 2) ~= n && size(O{i}, 1) == n
            O{i} = O{i}';
        end
        % Check that obstacle dimensions now match n.
        if size(O{i}, 2) ~= n
            error('Obstacle dimensions do not match the ellipsoid dimensionality.');
        end
        % Transpose obstacle to get an n x mₒ matrix, subtract bE, apply transformation, then transpose back.
        temp = (inv(E.DE) * E.AE') * (O{i}' - repmat(E.bE, 1, size(O{i},1)));
        O_transformed{i} = temp';  % Back to mₒ x n.
    end
end

% Use CasADi to compute the optimal halfspace direction for a given obstacle.
function ai = maximize_halfspace(Q, O)
    % Q is n x m and O is mₒ x n.
    n = size(Q, 1);
    
    % Define symbolic variable for halfspace direction using casadi.MX.sym.
    ai_sym = casadi.MX.sym('ai', n);
    
    % Define an objective function.
    % For demonstration, we maximize the negative sum of dot products with seed points.
    objective = -sum(dot(Q, repmat(ai_sym,1,size(Q,2))));
    
    % Define constraints: for each point in O, a'*point <= 0.
    constr = [];
    for j = 1:size(O, 1)
        constr = [constr; O(j,:) * ai_sym];
    end
    
    % Constraint bounds: each must be <= 0.
    lb_constr = -inf(size(constr));
    ub_constr = zeros(size(constr));
    
    % Set up options using nested structure fields.
    opts = struct();
    opts.ipopt = struct();
    opts.ipopt.print_level = 0;
    opts.print_time = 0;
    
    % Setup and solve the NLP with CasADi's nlpsol.
    nlp = struct('f', objective, 'x', ai_sym, 'g', constr);
    solver = casadi.nlpsol('solver', 'ipopt', nlp, opts);
    
    % Solve the optimization.
    sol = solver('x0', zeros(n,1), 'lbx', -inf(n,1), 'ubx', inf(n,1), ...
                 'lbg', lb_constr, 'ubg', ub_constr);
    ai = full(sol.x);
end

% Add a halfspace (represented by its normal vector) to the polytope P.
function P = add_halfspace(P, ai)
    % Placeholder: Append the halfspace vector (transposed for row format).
    % In a full implementation, P would store inequality constraints.
    P = [P; ai'];
end

% Update the polytope Pk by transforming the halfspaces back to the original space.
function Pk = update_polytope(P, E)
    % Placeholder: Recover the polytope using the ellipsoid's affine map.
    % Here we simply apply the inverse transformation.
    Pk = E.AE * P' + repmat(E.bE, 1, size(P,1));
end

% Compute the Maximum Volume Inscribed Ellipsoid (MVIE) for polytope P.
function Ek = compute_mvie(P)
    % Placeholder implementation.
    % A complete implementation would solve an optimization to find the largest ellipsoid
    % contained in the polytope defined by P.
    n = size(P, 2);  % Assuming P has columns equal to space dimension.
    Ek.AE = eye(n);
    Ek.DE = diag(rand(n,1) + 0.5);
    Ek.bE = zeros(n,1);
end

% Compute the "volume" of an ellipsoid (proportional to the product of semi-axis lengths).
function vol = volume(E)
    vol = prod(diag(E.DE));  % Simplified computation.
end

%% Plotting Function (for 2D)
function plot_firi_results(E_initial, Q, O, Pk)
    figure;
    hold on; grid on; axis equal;
    xlabel('X'); ylabel('Y');
    title('FIRI: Obstacles, Initial Seed Ellipsoid, and Final Polytope');
    
    % Plot obstacles.
    for i = 1:length(O)
        % O{i} is an m x 2 matrix (each row a point).
        obs = O{i};
        % Compute convex hull for obstacle boundary.
        k_obs = convhull(obs(:,1), obs(:,2));
        fill(obs(k_obs,1), obs(k_obs,2), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
    end
    
    % Plot seed points.
    scatter(Q(1,:), Q(2,:), 50, 'g', 'filled');
    
    % Plot initial ellipsoid (seed ellipsoid).
    plot_ellipse(E_initial, 'b');
    
    % Plot final polytope (assumed to be a set of vertices).
    % Pk is expected as a 2 x N matrix of vertices.
    Pk = Pk';
    k_poly = convhull(Pk(:,1), Pk(:,2));
    plot(Pk(k_poly,1), Pk(k_poly,2), 'k-', 'LineWidth', 2);
    
    legend('Obstacles','Seed Points','Initial Ellipsoid','Final Polytope');
    hold off;
end

% Utility function to plot a 2D ellipsoid.
function plot_ellipse(E, colorStr)
    % E is an ellipsoid structure with fields AE, DE, and bE.
    % For a 2D ellipsoid, we use the parametric equation:
    %   x = bE + AE * DE * [cos(t); sin(t)]
    t = linspace(0, 2*pi, 100);
    circle = [cos(t); sin(t)];  % 2 x 100
    ellipse = E.bE + E.AE * E.DE * circle;
    plot(ellipse(1,:), ellipse(2,:), colorStr, 'LineWidth', 2);
end


% Example inputs:
A = [1, 0; -1, 0];
b = [30; 20];
lb = [0, 0];
ub = [50, 50];
Q = rand(2, 10);              % 2 x 10 seed matrix (2D points)
O = {rand(5, 2), rand(4, 2)};   % Two obstacles; each is m x 2 (each row is a point)

% Run the algorithm.
Pk = firi_ell(A, b, lb, ub, Q, O);
