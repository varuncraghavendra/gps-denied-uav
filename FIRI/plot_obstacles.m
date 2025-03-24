function plot_nonoverlapping_convex_obstacles()
    % Define environment parameters
    env_size = [50, 50]; % Environment size (50m x 50m)
    side_length = 6; % Convex hull boundary (6m x 6m around seed)
    num_convex_obstacles = 10; % Number of convex obstacles
    min_spacing = 3; % Minimum spacing between obstacles

    % Generate a single random seed point
    seed = generate_random_seed(env_size);

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
    title('Non-Overlapping Convex Obstacles in 2D Space');

    % Plot convex obstacles
    for i = 1:length(convex_obstacles)
        k = convhull(convex_obstacles{i}(:,1), convex_obstacles{i}(:,2));
        fill(convex_obstacles{i}(k,1), convex_obstacles{i}(k,2), 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end

    % Draw the 6m x 6m boundary around the seed
    boundary_x = [seed(1) - side_length/2, seed(1) + side_length/2, ...
                  seed(1) + side_length/2, seed(1) - side_length/2, seed(1) - side_length/2];
    boundary_y = [seed(2) - side_length/2, seed(2) - side_length/2, ...
                  seed(2) + side_length/2, seed(2) + side_length/2, seed(2) - side_length/2];
    plot(boundary_x, boundary_y, 'b', 'LineWidth', 2);

    % Highlight the seed point
    plot(seed(1), seed(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

    hold off;
end

function seed = generate_random_seed(env_size)
    % Generate a random seed point in the environment
    seed = [rand * env_size(1), rand * env_size(2)];
end

function convex_obstacles = generate_nonoverlapping_obstacles(env_size, num_obstacles, min_spacing)
    % Generate non-overlapping convex obstacles using convex hulls
    convex_obstacles = cell(num_obstacles, 1);
    placed_centers = []; % Store placed obstacle centers
    
    for i = 1:num_obstacles
        num_points = randi([3, 8]); % Each obstacle has 3-8 points
        is_valid = false;
        
        % Ensure obstacles do not overlap and are spaced non-uniformly
        while ~is_valid
            center = [rand * env_size(1), rand * env_size(2)];
            
            % Check minimum spacing from other obstacles
            if isempty(placed_centers) || all(vecnorm(placed_centers - center, 2, 2) > min_spacing)
                is_valid = true;
                placed_centers = [placed_centers; center]; % Store the new center
                
                % Generate points around the center
                angles = linspace(0, 2*pi, num_points+1)' + rand * pi/4; % Random rotation
                angles(end) = []; % Remove last duplicate point
                radii = rand(num_points, 1) * 2 + 1; % Random radius between 1m-3m
                
                x = center(1) + radii .* cos(angles);
                y = center(2) + radii .* sin(angles);
                
                convex_obstacles{i} = [x, y];
            end
        end
    end
end

% Run the function
plot_nonoverlapping_convex_obstacles();
