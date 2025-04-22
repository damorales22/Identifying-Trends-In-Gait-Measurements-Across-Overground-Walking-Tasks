%% 3D Cube
close all
clear all

% Define the number of points in each dimension
n_points = 3;

% Generate the grid points for the cube
[x, y, z] = ndgrid(linspace(0, 1, n_points));

% Flatten the grids to create a list of points
x = x(:);
y = y(:);
z = z(:);

% Create the figure
figure;
hold on;

% Plot the points with increased size for thickness
scatter3(x, y, z, 200, 'filled', 'r'); % Increased size to 200 for thicker points

% Connect points to form a 3D grid (edges of the cube)
for i = 1:n_points
    for j = 1:n_points
        for k = 1:n_points
            % Connect along x-direction
            if i < n_points
                plot3([x((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       x(i*n_points^2 + (j-1)*n_points + k)], ...
                      [y((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       y(i*n_points^2 + (j-1)*n_points + k)], ...
                      [z((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       z(i*n_points^2 + (j-1)*n_points + k)], 'k', 'LineWidth', 1.5);
            end
            % Connect along y-direction
            if j < n_points
                plot3([x((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       x((i-1)*n_points^2 + j*n_points + k)], ...
                      [y((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       y((i-1)*n_points^2 + j*n_points + k)], ...
                      [z((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       z((i-1)*n_points^2 + j*n_points + k)], 'k', 'LineWidth', 1.5);
            end
            % Connect along z-direction
            if k < n_points
                plot3([x((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       x((i-1)*n_points^2 + (j-1)*n_points + (k+1))], ...
                      [y((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       y((i-1)*n_points^2 + (j-1)*n_points + (k+1))], ...
                      [z((i-1)*n_points^2 + (j-1)*n_points + k), ...
                       z((i-1)*n_points^2 + (j-1)*n_points + (k+1))], 'k', 'LineWidth', 1.5);
            end
        end
    end
end

% Custom solid arrows with thinner heads and no edges
arrow_length = 1.3; % Arrow length
arrow_head_size = 0.08; % Thinner size for the arrowhead
arrowhead_position = -0.05; % Variable to control the position of the arrowhead (0 for tip, > 0 for moving it away)

% X-axis (teal)
plot3([0, arrow_length], [0, 0], [0, 0], 'color', [0.25, 0.8, 0.8], 'LineWidth', 3); % Shaft
fill3([arrow_length-arrowhead_position, arrow_length-arrowhead_position-arrow_head_size, ...
       arrow_length-arrowhead_position-arrow_head_size], ...
      [0, arrow_head_size, -arrow_head_size], ...
      [0, 0, 0], [0.25, 0.8, 0.8], 'EdgeColor', 'none'); % Solid head with no edge

% Y-axis (lime green)
plot3([0, 0], [0, arrow_length], [0, 0], 'color', [0.6, 0.8, 0.6], 'LineWidth', 3); % Shaft
fill3([0, arrow_head_size, -arrow_head_size], ...
      [arrow_length-arrowhead_position, arrow_length-arrowhead_position-arrow_head_size, ...
       arrow_length-arrowhead_position-arrow_head_size], ...
      [0, 0, 0], [0.6, 0.8, 0.6], 'EdgeColor', 'none'); % Solid head with no edge

% Z-axis (royal blue)
plot3([0, 0], [0, 0], [0, arrow_length], 'color', [0.25, 0.4, 0.8], 'LineWidth', 3); % Shaft
fill3([0, 0, 0], ...
      [0, arrow_head_size, -arrow_head_size], ...
      [arrow_length-arrowhead_position, arrow_length-arrowhead_position-arrow_head_size, ...
       arrow_length-arrowhead_position-arrow_head_size], [0.25, 0.4, 0.8], 'EdgeColor', 'none'); % Solid head with no edge

% Set a 3D view angle
view(3);

% Remove axes and grid
axis equal;
axis off;
hold off;
