function plot_nonoverlapping_convex_obstacles_with_constraints(A, b, C, d, lb, ub)
    import iris.thirdParty.polytopes.*;
    
    % Define environment parameters
    env_size = [50, 50]; % Environment size (50m x 50m)
    num_convex_obstacles = 10; % Number of convex obstacles
    min_spacing = 3; % Minimum spacing between obstacles
    
    % Generate non-overlapping convex obstacles
    convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_convex_obstacles, min_spacing);
    
    % Plot the environment
    figure;
    hold on;
    grid on;
    axis equal;
    xlim([0, env_size(1)]);
    ylim([0, env_size(2)]);
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Non-Overlapping Convex Obstacles with Constraints');
    
    % Draw obstacles
    for i = 1:length(convex_obstacles)
        obs = convex_obstacles{i}';
        if size(obs, 2) > 2
            k = convhull(obs(1,:), obs(2,:));
        else
            k = [1,2,1];
        end
        patch(obs(1,k), obs(2,k), 'k', 'FaceColor', [.6,.6,.6], 'LineWidth', 0.1);
        plot(obs(1,k), obs(2,k), 'k', 'LineWidth', 2);
    end
    
    % Draw constraint boundaries if A and b are provided
    if ~isempty(A) && ~isempty(b)
        for j = 1:size(A,1)
            ai = A(j,:);
            bi = b(j);
            if ai(2) == 0
                x0 = [bi/ai(1); 0];
            else
                x0 = [0; bi/ai(2)];
            end
            u = [0,-1;1,0] * ai';
            pts = [x0 - 1000*u, x0 + 1000*u];
            plot(pts(1,:), pts(2,:), 'm--', 'LineWidth', 1.5);
        end
    end
    
    % Draw feasible region if A and b are provided
    if ~isempty(A) && ~isempty(b)
        V = lcon2vert(A, b);
        if ~isempty(V)
            k = convhull(V(:,1), V(:,2));
            plot(V(k,1), V(k,2), 'ro-', 'LineWidth', 2);
        end
    end
    
    % Draw ellipsoid constraint if C and d are provided
    if ~isempty(C) && ~isempty(d)
        th = linspace(0,2*pi,100);
        y = [cos(th);sin(th)];
        x = bsxfun(@plus, C*y, d);
        plot(x(1,:), x(2,:), 'b-', 'LineWidth', 2);
    end
    
    % Draw bounding box if lb and ub are provided
    if ~isempty(lb) && ~isempty(ub)
        plot([lb(1),ub(1),ub(1),lb(1),lb(1)], [lb(2),lb(2),ub(2),ub(2),lb(2)], 'k-');
        pad = (ub - lb) * 0.05;
        xlim([lb(1)-pad(1),ub(1)+pad(1)]);
        ylim([lb(2)-pad(2),ub(2)+pad(2)]);
    end
    
    axis off;
end

function convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_obstacles, min_spacing)
    convex_obstacles = cell(num_obstacles, 1);
    placed_centers = [];
    for i = 1:num_obstacles
        num_points = randi([3, 8]);
        is_valid = false;
        while ~is_valid
            center = [rand * env_size(1), rand * env_size(2)];
            if isempty(placed_centers) || all(vecnorm(placed_centers - center, 2, 2) > min_spacing)
                is_valid = true;
                placed_centers = [placed_centers; center];
                angles = linspace(0, 2*pi, num_points+1)' + rand * pi/4;
                angles(end) = [];
                radii = rand(num_points, 1) * 2 + 1;
                x = center(1) + radii .* cos(angles);
                y = center(2) + radii .* sin(angles);
                convex_obstacles{i} = [x, y];
            end
        end
    end
end
