%% Step Error Function
% this function inputs the marker locations (points struct), step times 
% (steps struct) and the M5 coordinates (trial struct and markerdata 
% struct) and outputs the matches struct containing Left and Right fields 
% each with (marker, timestamp, distance)

function [matches] = StepError(points, steps, markerdata, trialname, filepath) %points: red tape coordinates, steps: when steps are taken, markerdata: data from each sensor
    %% Inputs
    video = 0; % 0 if you don't want any recording, 1 if you want
        trial_for_video = 26; %If you are interested in record a video, specify the trial
        max_iteration = 15;   %Max of steps plotted before it clears the previous steps. Set a huge number to not take this into account
        pause_t = 0;          %Pause after plotting each step pair

    plot_longitudinal_error = 1;  %Alternate steps, longitudinal error. Choose this by DEFAULT as it is the only one working with video recording    
    plot_total_error = 0;         %Plot distance RIGHT step - Red tape (TOTAL error)
    plot_longitudinal_error2 = 0; %Plot distance RIGHT step - Infinite line (LONGITUDINAL error)
    
    %% Get Step Locations
    left = zeros(length(steps.left.frames),3);
    right = zeros(length(steps.right.frames),3);

    % left
    for i = 1:length(left(:,1))
        timestamp = steps.left.frames(i);
        step = markerdata.LMTH5(:, 1:3, timestamp);
        left(i, 1:3) = step; %position when each step was made
    end 
    % right
    for i = 1:length(right(:,1))
        timestamp = steps.right.frames(i);
        step = markerdata.RMTH5(:, 1:3, timestamp);
        right(i, 1:3) = step(1, :);
    end

    num = str2double(extractAfter(trialname, "Trial"));
    
    %% Identify straight lines and curves
    floor_markers_names = fieldnames(points);
    straight_points = [];
    curve_points = [];
    
    max_alignment_error = 50; %Max error (units) in the x axis we admit to clasify it as a straight horizontal line
    for i=1:length(fieldnames(points))
        if abs(abs(points.("FL"+i)(2,1)) - abs(points.("FL"+i)(2,2)))<max_alignment_error  %points.FLxx => first column: coordinates inner sensor; second column:coordinates external sensor
            straight_points{end+1} = floor_markers_names{i};
        else
            curve_points{end+1} = floor_markers_names{i};
        end
    end  
    
    %% Find Marker Match + Orthogonal Distance
    match_r = zeros(length(right(:,1)),3);
    match_l = zeros(length(left(:,1)),3);
    numFLmarkers = length(fieldnames(points)); %num of sensors used

    % function for dist between a point and a line (adapted to 2D; absolute value)
    function d = point_to_line(pt, v1, v2)
        %a = v1 - v2;
        %b = pt - v2;
        %d = norm(cross(a,b)) / norm(a);
        a = v2 - v1;  
        b = pt - v1;  
        d = abs(a(1) * b(2) - a(2) * b(1)) / norm(a);
    end

    function d = point_to_line_2(pt, v1, v2) %Not absolute: positive and negative values
        a = v2 - v1;  
        b = pt - v1;  
        d = (a(1) * b(2) - a(2) * b(1)) / norm(a);
    end

    % Get match_r for RIGHT foot
    %{
    for i=1:length(right(:,1))
        shortest = norm(points.FL1(:,3)' - right(i,:)); %For the first sensor. Column 3 represents the middle point between FL1a and b
        ind = 1;
        for j=2:numFLmarkers
            temp = norm(points.(['FL' num2str(j)])(:,3)' - right(i,:)); %From the second sensor
            if temp < shortest
                shortest = temp;
                ind = j;
            end
        end
        match_r(i,1) = steps.right.frames(i); %frame number of the step
        match_r(i,2) = ind; %index of the closest marker
        match_r(i,3) = point_to_line(right(i,:)', points.(['FL' num2str(ind)])(:,1), points.(['FL' num2str(ind)])(:,2)); %distance from the foot position to the closest marker
    end
    %}
    count_straight = 1;
    count_curve = 1;
    
    for i=1:length(right(:,1)) %for each step
        refpoint = right(i,1:2); 
        for j=1:numFLmarkers
            middle_tape = points.("FL" + num2str(j))(1:2,3)'; %get the x,y coordinates of the middle point of marker FLj
            distances(j,1) = sqrt((middle_tape(1,1) - refpoint(1,1)).^2 + (middle_tape(1,2) - refpoint(1,2)).^2); %get the distance between where the step takes place and the middle of every tape (FL segment)
        end
        [min_distance, idx] = min(distances); %we only care about idx (index) which represents the segment from which the step point is closest to. Then, we will use this to calculate the actual distance to the line between FlXa and FLXb.
        
        match_r(i,1) = steps.right.frames(i); %frame number of the step
        match_r(i,2) = idx; %index of the closest marker
        match_r(i,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest INFINITE line
        match_r(i,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2));

        if ismember("FL" + num2str(idx), straight_points) %Check if FLx is included in the straight section points
            match_r_straight(count_straight,1) = steps.right.frames(i); %frame number of the step
            match_r_straight(count_straight,2) = idx; %index of the closest marker
            match_r_straight(count_straight,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest INFINITE line
            match_r_straight(count_straight,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %not normalized error
            count_straight = count_straight + 1;
        elseif ismember("FL" + num2str(idx), curve_points)
            match_r_curve(count_curve,1) = steps.right.frames(i); %frame number of the step
            match_r_curve(count_curve,2) = idx; %index of the closest marker
            match_r_curve(count_curve,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest INFINITE line
            match_r_curve(count_curve,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2));
            count_curve = count_curve + 1;
        else
            error("FL" + num2str(idx) + " not included in any straight or curve section")
        end 
    end

    % Get matches for LEFT foot
    %{
    for i=1:length(left(:,1))
        shortest = norm(points.FL1(:,3)' - left(i,:));
        ind = 1;
        for j=2:numFLmarkers
            temp = norm(points.(['FL' num2str(j)])(:,3)' - left(i,:));
            if temp < shortest
                shortest = temp;
                ind = j;
            end
        end
        match_l(i,1) = steps.left.frames(i); %frame number of the step
        match_l(i,2) = ind; %index of the closest marker
        match_l(i,3) = point_to_line(left(i,:)', points.(['FL' num2str(ind)])(:,1), points.(['FL' num2str(ind)])(:,2)); %distance from the foot position to the closest marker
    end
    %}
    count_straight = 1;
    count_curve = 1;
    
    for i=1:length(left(:,1)) %for each step
        refpoint = left(i,1:2); 
        for j=1:numFLmarkers
            middle_tape = points.("FL" + num2str(j))(1:2,3)'; %get the x,y coordinates of the middle point of marker FLj
            distances(j,1) = sqrt((middle_tape(1,1) - refpoint(1,1)).^2 + (middle_tape(1,2) - refpoint(1,2)).^2); %get the distance between where the step takes place and the middle of every tape (FL segment)
        end
        [min_distance, idx] = min(distances); %we only care about idx (index) which represents the segment from which the step point is closest to. Then, we will use this to calculate the actual distance to the line between FlXa and FLXb.
        
        match_l(i,1) = steps.left.frames(i); %frame number of the step
        match_l(i,2) = idx; %index of the closest marker
        match_l(i,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest line
        match_l(i,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2));

        if ismember("FL" + num2str(idx), straight_points) %Check if FLx is included in the straight section points
            match_l_straight(count_straight,1) = steps.left.frames(i); %frame number of the step
            match_l_straight(count_straight,2) = idx; %index of the closest marker
            match_l_straight(count_straight,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest INFINITE line
            match_l_straight(count_straight,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %positive and negative errors
            count_straight = count_straight + 1;
        elseif ismember("FL" + num2str(idx), curve_points)
            match_l_curve(count_curve,1) = steps.left.frames(i); %frame number of the step
            match_l_curve(count_curve,2) = idx; %index of the closest marker
            match_l_curve(count_curve,3) = point_to_line(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2)); %distance from the foot position to the closest INFINITE line
            match_l_curve(count_curve,4) = point_to_line_2(refpoint', points.(['FL' num2str(idx)])(1:2,1), points.(['FL' num2str(idx)])(1:2,2));
            count_curve = count_curve + 1;
        else
            error("FL" + num2str(idx) + " not included in any straight or curve section")
        end
    end    
    clear idx;
    
    %% Define data
    
    % means_r = zeros(1, numFLmarkers);
    % stds_r = zeros(1, numFLmarkers);
    % means_l = zeros(1, numFLmarkers);
    % stds_l = zeros(1, numFLmarkers);

    % for i = 1:numFLmarkers  % calculate and analyze the error for each marker separately and then combine the results
    %     isolated = match_r(match_r(:, 2) == i, :); % isolate matches for the ith marker
    % 
    %     mean_value = mean(abs(isolated(:, 3))); % calculate mean for marker
    %     std_value = std(abs(isolated(:, 3))); % calculate standard deviation
    % 
    %     % Store unique value and mean in 'means' array
    %     means_r(i) = mean_value; %mean orthogonal distance for each marker
    %     stds_r(i) = std_value;   %std of the orthogonal distance for each marker
    % 
    %     isolated = match_l(match_l(:, 2) == i, :); % isolate matches for the ith marker
    % 
    %     mean_value = mean(abs(isolated(:, 3))); % calculate mean for marker
    %     std_value = std(abs(isolated(:, 3))); % calculate standard deviation
    % 
    %     % Store unique value and mean in 'means' array
    %     means_l(i) = mean_value;
    %     stds_l(i) = std_value; 
    % end

    % % For all sections - Absolute error
    % overall_mean_l = mean(means_l(~isnan(means_l)));
    % overall_mean_r = mean(means_r(~isnan(means_r)));
    % overall_mean = (overall_mean_r + overall_mean_l)/2; %overall mean of both feet combined
    % overall_std_l = std(means_l(~isnan(means_l)));
    % overall_std_r = std(means_r(~isnan(means_l)));
    % overall_std = (overall_std_r + overall_std_l)/2;
    
    match_l(any(isnan(match_l), 2), :) = [];
    match_r(any(isnan(match_r), 2), :) = [];
    match_l_straight(any(isnan(match_l_straight), 2), :) = [];
    match_r_straight(any(isnan(match_r_straight), 2), :) = [];

    % For all sections - Absolute error
    overall_mean_l = mean(match_l(:,3));
    overall_mean_r = mean(match_r(:,3));
    overall_mean = (overall_mean_r + overall_mean_l)/2; %overall mean of both feet combined
    overall_std_l = std(match_l(:,3));
    overall_std_r = std(match_r(:,3));
    overall_std = (overall_std_r + overall_std_l)/2;

    % For all sections - Signed error
    overall_mean_l_signed = mean(match_l(:,4));
    overall_mean_r_signed = mean(match_r(:,4));
    overall_mean_signed = (overall_mean_r_signed + overall_mean_l_signed)/2; %overall mean of both feet combined
    overall_std_l_signed = std(match_l(:,4));
    overall_std_r_signed = std(match_r(:,4));
    overall_std_signed = (overall_std_r_signed + overall_std_l_signed)/2;

    % Straight sections - Absolute error
    overall_mean_l_straight = mean(match_l_straight(:,3));
    overall_mean_r_straight = mean(match_r_straight(:,3));
    overall_mean_straight = (overall_mean_r_straight + overall_mean_l_straight)/2; %overall mean of both feet combined
    overall_median_straight = median([match_l_straight(:,3); match_r_straight(:,3)]);
    overall_std_l_straight = std(match_l_straight(:,3));
    overall_std_r_straight = std(match_r_straight(:,3));
    overall_std_straight = (overall_std_r_straight + overall_std_l_straight)/2;

    % Straight sections - Signed error
    overall_mean_l_straight_signed = mean(match_l_straight(:,4));
    overall_mean_r_straight_signed = mean(match_r_straight(:,4));
    overall_mean_straight_signed = (overall_mean_r_straight_signed + overall_mean_l_straight_signed)/2; %overall mean of both feet combined
    overall_std_l_straight_signed = std(match_l_straight(:,4));
    overall_std_r_straight_signed = std(match_r_straight(:,4));
    overall_std_straight_signed = (overall_std_r_straight_signed + overall_std_l_straight_signed)/2;

    % Curves - Absolute error
    overall_mean_l_curve = mean(match_l_curve(:,3));
    overall_mean_r_curve = mean(match_r_curve(:,3));
    overall_mean_curve = (overall_mean_r_curve + overall_mean_l_curve)/2; %overall mean of both feet combined
    overall_std_l_curve = std(match_l_curve(:,3));
    overall_std_r_curve = std(match_r_curve(:,3));
    overall_std_curve = (overall_std_r_curve + overall_std_l_curve)/2;

    % Curves - Signed error
    overall_mean_l_curve_signed = mean(match_l_curve(:,4));
    overall_mean_r_curve_signed = mean(match_r_curve(:,4));
    overall_mean_curve_signed = (overall_mean_r_curve_signed + overall_mean_l_curve_signed)/2; %overall mean of both feet combined
    overall_std_l_curve_signed = std(match_l_curve(:,4));
    overall_std_r_curve_signed = std(match_r_curve(:,4));
    overall_std_curve_signed = (overall_std_r_curve_signed + overall_std_l_curve_signed)/2;

%% Store data

    matches = struct();

    matches.right.match = match_r;
    matches.right.match_straight = match_r_straight;
    matches.right.match_curve = match_r_curve;

    %matches.right.means = means_r;
    matches.right.overall_mean = overall_mean_r;
    matches.right.overall_std = overall_std_r;
    matches.right.straight_mean = overall_mean_r_straight;
    matches.right.straight_std = overall_std_r_straight;
    matches.right.curves_mean = overall_mean_r_curve;
    matches.right.curves_std = overall_std_r_curve;

    matches.left.match = match_l;
    matches.left.match_straight = match_l_straight;
    matches.left.match_curve = match_l_curve;
    
    %matches.left.means = means_l;
    matches.left.overall_mean = overall_mean_l;
    matches.left.overall_std = overall_std_l;
    matches.left.straight_mean = overall_mean_l_straight;
    matches.left.straight_std = overall_std_l_straight;
    matches.left.curves_mean = overall_mean_l_curve;
    matches.left.curves_std = overall_std_l_curve;

    matches.overall.overall_mean_absolute = overall_mean;
    matches.overall.overall_std_absolute = overall_std;
    matches.overall.overall_mean_signed = overall_mean_signed;
    matches.overall.overall_std_signed = overall_std_signed;

    matches.overall.straight_mean_absolute = overall_mean_straight;
    matches.overall.straight_std_absolute = overall_std_straight;
    matches.overall.straight_median = overall_median_straight;
    matches.overall.curve_mean_absolute = overall_mean_curve;
    matches.overall.curve_std_absolute = overall_std_curve;

    matches.overall.straight_mean_signed = overall_mean_straight_signed;
    matches.overall.straight_std_signed = overall_std_straight_signed;
    matches.overall.curve_mean_signed = overall_mean_curve_signed;
    matches.overall.curve_std_signed = overall_std_curve_signed;

    matches.straight.markers = straight_points;
    matches.curve.markers = curve_points;

%% Plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT DISTANCE STEP - RED TAPE (TOTAL ERROR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot_total_error == 1
        randomColor = rand(1, 3);  % Random RGB values between 0 and 1
        clf; %To clear what it was plotted in the previous loop
        hold on
    
        % Plot the floor markers (points and lines between them)
        for j = 1:length(fieldnames(points))
            plot(points.("FL" + num2str(j))(1, 1:2), points.("FL" + num2str(j))(2, 1:2), 'o'); %Plot the points of the floor markers
    
            % Plot a line between the two points of the floor marker
            line_x = points.("FL" + num2str(j))(1, 1:2);
            line_y = points.("FL" + num2str(j))(2, 1:2);
            plot(line_x, line_y, 'Color', 'k');  % Plot the line in black
        end
    
        %legend()
    
        % Plot the steps and the red line from the step to the line between markers
        for j = 1:length(matches.right.match)
            timestamp = matches.right.match(j, 1);
            step = markerdata.RMTH5(:, 1:2, timestamp);  % Step point
            plot(step(1, 1), step(1, 2), 'o', 'MarkerFaceColor', randomColor, 'MarkerEdgeColor', randomColor);  % Plot the step point
    
            % Find the closest floor marker segment
            idx = matches.right.match(j, 2);  % Index of the closest floor marker
    
            % Get the endpoints of the closest floor marker line segment
            P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
            P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
    
            % Inline calculation to project the step point onto the line defined by P1 and P2
            v = P2 - P1;  % Vector from P1 to P2
            w = step(1, 1:2) - P1';  % Vector from P1 to the step point
            t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line defined by P1 and P2
            t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
            proj_point = P1 + t * v; % Calculate the projection point
    
            %Plot a red line from the step point to the projection on the segment
            plot([step(1, 1), proj_point(1)], [step(1, 2), proj_point(2)], 'r-', 'LineWidth', 2); 
    
            %Calculate the Euclidean distance between the step point and the projection point (closest point of the red tape)
            distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
            matches.right.match(j, 5) = distance;
        end
    
        axis equal
        axis square;
        grid on;
        title('Right foot and floor markers distribution')
        xlim([-6500, 3500]);
        ylim([-4000, 6000]);
        hold off
    
        % Initialize distances for the left steps
        for j = 1:length(matches.left.match)
            timestamp = matches.left.match(j, 1);
            step = markerdata.LMTH5(:, 1:2, timestamp);  % Step point
    
            % Find the closest floor marker segment
            idx = matches.left.match(j, 2);  % Index of the closest floor marker
    
            % Get the endpoints of the closest floor marker line segment
            P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
            P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
    
            % Inline calculation to project the step point onto the line defined by P1 and P2
            v = P2 - P1;  % Vector from P1 to P2
            w = step(1, 1:2) - P1';  % Vector from P1 to the step point
            t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line defined by P1 and P2
    
            % No need to clamp t since we're calculating projection onto the line
            proj_point = P1 + t * v; % Calculate the projection point
    
            % Calculate the Euclidean distance between the step point and the projection point
            distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
    
            % Store the distance in the match structure
            matches.left.match(j, 5) = distance;  % Optional: store it in a different column 
        end
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT DISTANCE STEP - INFINITE LINE (LONGITUDINAL ERROR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot_longitudinal_error2 == 1
        clf; %To clear what it was plotted in the previous loop
        hold on;
        axis equal;
        axis square;
        %title('Distance to Infinite Line Defined by Floor Markers - Right Foot');
        title('Distance to Infinite Line Defined by Floor Markers - ' + string(trialname) + ' - Right Foot');
        grid on;
        xlim([-6500, 3500]);
        ylim([-4000, 6000]);
    
        % Plot the floor markers (points and lines between them)
        for j = 1:length(fieldnames(points))
            % Plot the points of the floor markers
            plot(points.("FL" + num2str(j))(1, 1:2), points.("FL" + num2str(j))(2, 1:2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
    
            % Get the endpoints of the floor marker line segment
            P1 = points.("FL" + num2str(j))(1:2, 1);  % First point of the segment
            P2 = points.("FL" + num2str(j))(1:2, 2);  % Second point of the segment
    
            % Calculate the direction vector and extend the line
            v = P2 - P1;  % Vector from P1 to P2
    
            % Extend the line in both directions
            extended_line_x = [P1(1) - 0.75 * v(1), P2(1) + 0.75 * v(1)];  
            extended_line_y = [P1(2) - 0.75 * v(2), P2(2) + 0.75 * v(2)];  
    
            % Plot the infinite line in black dashed style
            plot(extended_line_x, extended_line_y, 'k--');  
    
            % Determine which point is exterior and add the floor marker label
            label_offset = 225;  % Increased this value for much more spacing
    
            % Choose the exterior point for label placement (using P2 here for exterior)
            label_position_right = P2 + [label_offset, 0];  % Place label to the right
            label_position_left = P2 + [-label_offset, 0];   % Place label to the left
    
            % Add floor marker label next to the exterior point (to the right)
            text(label_position_right(1), label_position_right(2), "FL" + num2str(j), ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
                'FontSize', 10, 'Color', 'k');
    
            % Optionally, if you want to place it to the left instead, uncomment the line below
            % text(label_position_left(1), label_position_left(2), "FL" + num2str(j), ...
            %     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
            %     'FontSize', 10, 'Color', 'k');
        end   
    
        % legend('Floor Markers', 'Location', 'best');
        randomColor = rand(1, 3);  % Random RGB values between 0 and 1
    
        count_straight = 1;
        count_curve = 1;
        % Independent plot for the distance to the infinite line
        for j = 1:length(matches.right.match)
            timestamp = matches.right.match(j, 1);
            step = markerdata.RMTH5(:, 1:2, timestamp);  % Step point
    
            % Find the closest floor marker segment
            idx = matches.right.match(j, 2);  % Index of the closest floor marker
    
            % Get the endpoints of the closest floor marker line segment
            P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
            P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
    
            % Calculate the projection onto the infinite line defined by P1 and P2
            v = P2 - P1;  % Vector from P1 to P2
            w = step(1, 1:2) - P1';  % Vector from P1 to the step point
            t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line (without clamping t)
            proj_point = P1 + t * v;  % Calculate the projection point on the infinite line
    
            % Plot the step point
            plot(step(1, 1), step(1, 2), 'o', 'MarkerFaceColor', randomColor, 'MarkerEdgeColor', randomColor);
    
            % Plot a red line from the step point to the projection on the infinite line
            plot([step(1, 1), proj_point(1)], [step(1, 2), proj_point(2)], 'r-', 'LineWidth', 2);
    
            % To save the distance from the match to the segment
            t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
            proj_point = P1 + t * v;
            distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
            matches.right.match(j, 5) = distance;
    
            if matches.right.match(j,2) == matches.right.match_straight(count_straight,2)
                matches.right.match_straight(count_straight, 5) = distance;
                if count_straight < length(matches.right.match_straight(:,3))
                    count_straight = count_straight + 1;
                end
            else
                matches.right.match_curve(count_curve, 5) = distance;
                count_curve = count_curve + 1;
            end   
    
        end %BREAK here to plot step by step
    
        hold off;
    
        %Comment until the next for if you don't want to plot the left foot (+ the two last plot lines in the next for)
        % figure
        % hold on
        % %To plot the floor lines and markers: comment if you don't want to plot the left foot steps
        % for j = 1:length(fieldnames(points))
        %     % Plot the points of the floor markers
        %     plot(points.("FL" + num2str(j))(1, 1:2), points.("FL" + num2str(j))(2, 1:2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
        % 
        %     % Get the endpoints of the floor marker line segment
        %     P1 = points.("FL" + num2str(j))(1:2, 1);  % First point of the segment
        %     P2 = points.("FL" + num2str(j))(1:2, 2);  % Second point of the segment
        % 
        %     % Calculate the direction vector and extend the line
        %     v = P2 - P1;  % Vector from P1 to P2
        % 
        %     % Extend the line in both directions
        %     extended_line_x = [P1(1) - 0.75 * v(1), P2(1) + 0.75 * v(1)];  
        %     extended_line_y = [P1(2) - 0.75 * v(2), P2(2) + 0.75 * v(2)];  
        % 
        %     % Plot the infinite line in black dashed style
        %     plot(extended_line_x, extended_line_y, 'k--');  
        % 
        %     % Determine which point is exterior and add the floor marker label
        %     label_offset = 225;  % Increased this value for much more spacing
        % 
        %     % Choose the exterior point for label placement (using P2 here for exterior)
        %     label_position_right = P2 + [label_offset, 0];  % Place label to the right
        %     label_position_left = P2 + [-label_offset, 0];   % Place label to the left
        % 
        %     % Add floor marker label next to the exterior point (to the right)
        %     text(label_position_right(1), label_position_right(2), "FL" + num2str(j), ...
        %         'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
        %         'FontSize', 10, 'Color', 'k');
        % end   
        % 
        % axis equal;
        % axis square;
        % title('Distance to Infinite Line Defined by Floor Markers - Left Foot');
        % grid on;
        % xlim([-6500, 3500]);
        % ylim([-4000, 6000]);
    
        % Initialize distances for the left steps (useful just if you wanna plot it too. Otherwise, it doesn't do anything)
        count_straight = 1;
        count_curve = 1;
        for j = 1:length(matches.left.match)
            timestamp = matches.left.match(j, 1);
            step = markerdata.LMTH5(:, 1:2, timestamp);  % Step point
    
            % Find the closest floor marker segment
            idx = matches.left.match(j, 2);  % Index of the closest floor marker
    
            % Get the endpoints of the closest floor marker line segment
            P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
            P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
    
            % Inline calculation to project the step point onto the line defined by P1 and P2
            v = P2 - P1;  % Vector from P1 to P2
            w = step(1, 1:2) - P1';  % Vector from P1 to the step point
            t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line defined by P1 and P2
    
            % No need to clamp t since we're calculating projection onto the line
            proj_point = P1 + t * v; % Calculate the projection point
    
            %To save the distance from the match to the segment
            %t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
            proj_point = P1 + t * v;
            distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
            matches.left.match(j, 5) = distance;
    
            if matches.left.match(j,2) == matches.left.match_straight(count_straight,2)
                matches.left.match_straight(count_straight, 5) = distance;
                if count_straight < length(matches.left.match_straight(:,3))
                    count_straight = count_straight + 1;
                end
            else
                matches.left.match_curve(count_curve, 5) = distance;
                count_curve = count_curve + 1;
            end 
    
            % % Plot the step point
            % plot(step(1, 1), step(1, 2), 'o', 'MarkerFaceColor', randomColor, 'MarkerEdgeColor', randomColor);
            % hold on;  % Keep previous plot data
            % 
            % % Plot a red line from the step point to the projection on the infinite line
            % plot([step(1, 1), proj_point(1)], [step(1, 2), proj_point(2)], 'r-', 'LineWidth', 2);
        end
        %hold off;
    end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALTERNATE RIGHT AND LEFT STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot_longitudinal_error == 1    
        if video == 0 || num ~= trial_for_video %To plot everything NOT recording the video
            clf; % Clear previous plots
            hold on;
            axis equal;
            axis square;
            title('Distance to Infinite Line Defined by Floor Markers - ' + string(trialname));
            
            %grid on;
            ax = gca;  % Get the current axes
            ax.GridColor = [0, 0, 0];  % Set grid color to black
            ax.GridAlpha = 0.3;  % Make the grid lines less transparent (darker)
    
            xticks(-650:50:350); 
            yticks(-400:50:600);
            axis equal; 
            xlim([-650, 350]); 
            ylim([-400, 600]);  
    
            % Plot the floor markers (points and lines between them)
            for j = 1:length(fieldnames(points))
                % Plot the points of the floor markers
                plot((points.("FL" + num2str(j))(1, 1:2))/10, (points.("FL" + num2str(j))(2, 1:2))/10, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
        
                % Get the endpoints of the floor marker line segment
                P1 = points.("FL" + num2str(j))(1:2, 1);  % First point of the segment
                P2 = points.("FL" + num2str(j))(1:2, 2);  % Second point of the segment
        
                % Calculate the direction vector and extend the line
                v = P2 - P1;  % Vector from P1 to P2
        
                % Extend the line in both directions
                extended_line_x = [P1(1) - 0.75 * v(1), P2(1) + 0.75 * v(1)];  
                extended_line_y = [P1(2) - 0.75 * v(2), P2(2) + 0.75 * v(2)];  
        
                % Plot the infinite line in black dashed style
                plot(extended_line_x/10, extended_line_y/10, 'k--');  
        
                % Add floor marker label
                label_offset = 225;  % Increased this value for more spacing
                label_position = P2 + [label_offset, 0];  % Place label to the right
                text(label_position(1)/10, label_position(2)/10, "FL" + num2str(j), ...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
                    'FontSize', 10, 'Color', 'k');
            end   
        
            % Initialize colors for plotting
            rightColor = 'r';   % Color for right steps (red)
            leftColor = 'b';    % Color for left steps (blue)
        
            % Initialize counters for right and left steps
            count_straight_right = 1;
            count_curve_right = 1;
            count_straight_left = 1;
            count_curve_left = 1;
        
            % Get the total number of steps for both sides
            numRightSteps = length(matches.right.match);
            numLeftSteps = length(matches.left.match);
            maxSteps = max(numRightSteps, numLeftSteps);
        
            % Alternate plotting of right and left steps
            for j = 1:maxSteps
                % Plot right steps if they exist
                if j <= numRightSteps
                    timestamp = matches.right.match(j, 1);
                    step = markerdata.RMTH5(:, 1:2, timestamp);  % Step point
        
                    % Find the closest floor marker segment
                    idx = matches.right.match(j, 2);  % Index of the closest floor marker
        
                    % Get the endpoints of the closest floor marker line segment
                    P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
                    P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
        
                    % Calculate the projection onto the infinite line defined by P1 and P2
                    v = P2 - P1;  % Vector from P1 to P2
                    w = step(1, 1:2) - P1';  % Vector from P1 to the step point
                    t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line
                    proj_point = P1 + t * v;  % Calculate the projection point on the infinite line
        
                    % Plot the right step point
                    plot(step(1, 1)/10, step(1, 2)/10, 'o', 'MarkerFaceColor', rightColor, 'MarkerEdgeColor', rightColor);
        
                    % Plot a red line from the step point to the projection on the infinite line
                    plot([step(1, 1)/10, proj_point(1)/10], [step(1, 2)/10, proj_point(2)/10], 'r-', 'LineWidth', 2);
        
                    % Calculate distance
                    t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
                    proj_point = P1 + t * v;
                    distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
                    matches.right.match(j, 5) = distance;
        
                    if matches.right.match(j, 2) == matches.right.match_straight(count_straight_right, 2)
                        matches.right.match_straight(count_straight_right, 5) = distance;
                        if count_straight_right < length(matches.right.match_straight(:, 3))
                            count_straight_right = count_straight_right + 1;
                        end
                    else
                        matches.right.match_curve(count_curve_right, 5) = distance;
                        count_curve_right = count_curve_right + 1;
                    end   
                end %BREAK HERE for alternating steps
        
                % Plot left steps if they exist
                if j <= numLeftSteps
                    timestamp = matches.left.match(j, 1);
                    step = markerdata.LMTH5(:, 1:2, timestamp);  % Step point
        
                    % Find the closest floor marker segment
                    idx = matches.left.match(j, 2);  % Index of the closest floor marker
        
                    % Get the endpoints of the closest floor marker line segment
                    P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
                    P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
        
                    % Calculate the projection onto the infinite line defined by P1 and P2
                    v = P2 - P1;  % Vector from P1 to P2
                    w = step(1, 1:2) - P1';  % Vector from P1 to the step point
                    t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line
                    proj_point = P1 + t * v;  % Calculate the projection point on the infinite line
        
                    % Plot the left step point
                    plot(step(1, 1)/10, step(1, 2)/10, 'o', 'MarkerFaceColor', leftColor, 'MarkerEdgeColor', leftColor); %BREAK HERE
        
                    % Plot a blue line from the step point to the projection on the infinite line
                    plot([step(1, 1)/10, proj_point(1)/10], [step(1, 2)/10, proj_point(2)/10], 'b-', 'LineWidth', 2);
        
                    % Calculate distance
                    t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
                    proj_point = P1 + t * v;
                    distance = sqrt((step(1, 1) - proj_point(1))^2 + (step(1, 2) - proj_point(2))^2); 
                    matches.left.match(j, 5) = distance;
        
                    if matches.left.match(j, 2) == matches.left.match_straight(count_straight_left, 2)
                        matches.left.match_straight(count_straight_left, 5) = distance;
                        if count_straight_left < length(matches.left.match_straight(:, 3))
                            count_straight_left = count_straight_left + 1;
                        end
                    else
                        matches.left.match_curve(count_curve_left, 5) = distance;
                        count_curve_left = count_curve_left + 1;
                    end   
                end 
            end %BREAK HERE
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VIDEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    
        if video == 1 && num ==trial_for_video
            clf;
            % Create a VideoWriter object
            videoFileName = filepath + 'Evolution - Trial00' + trial_for_video + '.avi';  % Specify the video file name
            videoWriter = VideoWriter(videoFileName); % Create a VideoWriter object
            videoWriter.FrameRate = 1; % Set the frame rate to 1 frame per second
            open(videoWriter); % Open the video file for writing
                  
            % Initialize colors for plotting
            rightColor = 'r';   % Color for right steps (red)
            leftColor = 'b';    % Color for left steps (blue)
            
            % Initialize counters for right and left steps
            count_straight_right = 1;
            count_curve_right = 1;
            count_straight_left = 1;
            count_curve_left = 1;
            
            % Buffer to keep track of previous iterations for clearing lines and dots
            previousRight = cell(1, 14); % Buffer for right step data
            previousLeft = cell(1, 14); % Buffer for left step data
            iterationBuffer = 1; % Current iteration counter
            
            % Get the total number of steps for both sides
            numRightSteps = length(matches.right.match);
            numLeftSteps = length(matches.left.match);
            maxSteps = max(numRightSteps, numLeftSteps);
            
            figure('units', 'normalized', 'outerposition', [0 0 1 1]); % Maximized figure window
            % Alternate plotting of right and left steps
            for j=1:maxSteps
                
                if iterationBuffer == 1
                    % Plot the floor markers (points and lines between them)                   
                    hold on;
                    axis equal;
                    axis square;
                    title('Distance to Infinite Line Defined by Floor Markers - ' + string(trialname));
    
                    grid on;
                    ax = gca;  % Get the current axes
                    ax.GridColor = [0, 0, 0];  % Set grid color to black
                    ax.GridAlpha = 0.2;  % Make the grid lines less transparent (darker)
    
                    xlabel("Distance (cm)")
                    ylabel("Distance (cm)")
                    xticks(-650:50:350); % X-axis ticks every 50 units
                    yticks(-400:50:600); % Y-axis ticks every 50 units
                    axis equal; % Or you can use axis square
                    xlim([-650, 350]);
                    ylim([-400, 600]);
    
                    for k = 1:length(fieldnames(points))
                        % Plot the points of the floor markers
                        plot((points.("FL" + num2str(k))(1, 1:2))/10, (points.("FL" + num2str(k))(2, 1:2))/10, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); % /10 to put it into cm
                    
                        % Get the endpoints of the floor marker line segment
                        P1 = points.("FL" + num2str(k))(1:2, 1);  % First point of the segment
                        P2 = points.("FL" + num2str(k))(1:2, 2);  % Second point of the segment
                    
                        % Calculate the direction vector and extend the line
                        v = P2 - P1;  % Vector from P1 to P2
                    
                        % Extend the line in both directions
                        extended_line_x = [P1(1) - 0.75 * v(1), P2(1) + 0.75 * v(1)];  
                        extended_line_y = [P1(2) - 0.75 * v(2), P2(2) + 0.75 * v(2)];  
                    
                        % Plot the infinite line in black dashed style
                        plot(extended_line_x/10, extended_line_y/10, 'k--');  
                    
                        % Add floor marker label
                        label_offset = 225;  % Increased this value for more spacing
                        label_position = P2 + [label_offset, 0];  % Place label to the right
                        text(label_position(1)/10, label_position(2)/10, "FL" + num2str(k), ...  % /10 to cm
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
                            'FontSize', 10, 'Color', 'k');
                    end 
                end
    
                % Plot right steps if they exist
                if j <= numRightSteps
                    timestamp = matches.right.match(j, 1);
                    step_r = markerdata.RMTH5(:, 1:2, timestamp);  % Step point
            
                    % Find the closest floor marker segment
                    idx = matches.right.match(j, 2);  % Index of the closest floor marker
            
                    % Get the endpoints of the closest floor marker line segment
                    P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
                    P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
            
                    % Calculate the projection onto the infinite line defined by P1 and P2
                    v = P2 - P1;  % Vector from P1 to P2
                    w = step_r(1, 1:2) - P1';  % Vector from P1 to the step point
                    t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line
                    proj_point_r = P1 + t * v;  % Calculate the projection point on the infinite line
                  
                    % % Plot the right step point
                    % plot(step_r(1, 1)/10, step_r(1, 2)/10, 'o', 'MarkerFaceColor', rightColor, 'MarkerEdgeColor', rightColor); % /10 to convert it to cm
                    % 
                    % % Plot a red line from the step point to the projection on the infinite line
                    % plot([step_r(1, 1)/10, proj_point_r(1)/10], [step_r(1, 2)/10, proj_point_r(2)/10], 'r-', 'LineWidth', 2);
            
                    % Calculate distance
                    t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
                    proj_point_r_2 = P1 + t * v;
                    distance = sqrt((step_r(1, 1) - proj_point_r_2(1))^2 + (step_r(1, 2) - proj_point_r_2(2))^2); 
                    matches.right.match(j, 5) = distance;
            
                    if matches.right.match(j, 2) == matches.right.match_straight(count_straight_right, 2)
                        matches.right.match_straight(count_straight_right, 5) = distance;
                        if count_straight_right < length(matches.right.match_straight(:, 3))
                            count_straight_right = count_straight_right + 1;
                        end
                    else
                        matches.right.match_curve(count_curve_right, 5) = distance;
                        count_curve_right = count_curve_right + 1;
                    end   
                end %BREAK HERE
            
                % frame = getframe(gcf); % Get the current figure frame
                % writeVideo(videoWriter, frame); % Write the frame to the video
                % pause(pause_t);
            
                % Plot left steps if they exist
                if j <= numLeftSteps
                    timestamp = matches.left.match(j, 1);
                    step_l = markerdata.LMTH5(:, 1:2, timestamp);  % Step point
            
                    % Find the closest floor marker segment
                    idx = matches.left.match(j, 2);  % Index of the closest floor marker
            
                    % Get the endpoints of the closest floor marker line segment
                    P1 = points.("FL" + num2str(idx))(1:2, 1);  % First point of the segment
                    P2 = points.("FL" + num2str(idx))(1:2, 2);  % Second point of the segment
            
                    % Calculate the projection onto the infinite line defined by P1 and P2
                    v = P2 - P1;  % Vector from P1 to P2
                    w = step_l(1, 1:2) - P1';  % Vector from P1 to the step point
                    t = dot(w, v) / dot(v, v);  % Project the step point onto the infinite line
                    proj_point_l = P1 + t * v;  % Calculate the projection point on the infinite line
            
                    % Check dimensions before storing in the buffer
                    if size(step_l, 2) == 2 && size(proj_point_l, 1) == 1 && size(proj_point_l, 2) == 2
                        if iterationBuffer <= 14
                            previousLeft{iterationBuffer} = [step_l(1, :), proj_point_l]; % Store step and projection
                        end
                    end
            
                    % % Plot the left step point
                    % plot(step_l(1, 1)/10, step_l(1, 2)/10, 'o', 'MarkerFaceColor', leftColor, 'MarkerEdgeColor', leftColor);
                    % 
                    % % Plot a blue line from the step point to the projection on the infinite line
                    % plot([step_l(1, 1)/10, proj_point_l(1)/10], [step_l(1, 2)/10, proj_point_l(2)/10], 'b-', 'LineWidth', 2); %to cm
            
                    % Calculate distance
                    t = max(0, min(1, t)); % Clamp t to ensure the projection is within the segment
                    proj_point_l_2 = P1 + t * v;
                    distance = sqrt((step_l(1, 1) - proj_point_l_2(1))^2 + (step_l(1, 2) - proj_point_l_2(2))^2); 
                    matches.left.match(j, 5) = distance;
            
                    if matches.left.match(j, 2) == matches.left.match_straight(count_straight_left, 2)
                        matches.left.match_straight(count_straight_left, 5) = distance;
                        if count_straight_left < length(matches.left.match_straight(:, 3))
                            count_straight_left = count_straight_left + 1;
                        end
                    else
                        matches.left.match_curve(count_curve_left, 5) = distance;
                        count_curve_left = count_curve_left + 1;
                    end
                end
    
                % Plot right and left steps, taking into account the first step side
                if steps.right.secs(1,1)<steps.left.secs(1,1) %the right step goes first    
                        first_step = "right"; %variable just for control
    
                        % Plot the right step point
                        plot(step_r(1, 1)/10, step_r(1, 2)/10, 'o', 'MarkerFaceColor', rightColor, 'MarkerEdgeColor', rightColor); % /10 to convert it to cm
                        % Plot a red line from the step point to the projection on the infinite line
                        plot([step_r(1, 1)/10, proj_point_r(1)/10], [step_r(1, 2)/10, proj_point_r(2)/10], 'r-', 'LineWidth', 2);
    
                        frame = getframe(gcf); % Get the current figure frame
                        writeVideo(videoWriter, frame); % Write the frame to the video
                        pause(pause_t);
                
                        % Plot the left step point
                        plot(step_l(1, 1)/10, step_l(1, 2)/10, 'o', 'MarkerFaceColor', leftColor, 'MarkerEdgeColor', leftColor);
                        % Plot a blue line from the step point to the projection on the infinite line
                        plot([step_l(1, 1)/10, proj_point_l(1)/10], [step_l(1, 2)/10, proj_point_l(2)/10], 'b-', 'LineWidth', 2); %to cm
    
                        frame = getframe(gcf); % Get the current figure frame
                        writeVideo(videoWriter, frame); % Write the frame to the video
                        pause(pause_t);
    
                    else %left step goes first
                        first_step = "left";  %variable just for control
    
                        % Plot the left step point
                        plot(step_l(1, 1)/10, step_l(1, 2)/10, 'o', 'MarkerFaceColor', leftColor, 'MarkerEdgeColor', leftColor);
                        % Plot a blue line from the step point to the projection on the infinite line
                        plot([step_l(1, 1)/10, proj_point_l(1)/10], [step_l(1, 2)/10, proj_point_l(2)/10], 'b-', 'LineWidth', 2); %to cm
    
                        frame = getframe(gcf); % Get the current figure frame
                        writeVideo(videoWriter, frame); % Write the frame to the video
                        pause(pause_t);
        
                        % Plot the right step point
                        plot(step_r(1, 1)/10, step_r(1, 2)/10, 'o', 'MarkerFaceColor', rightColor, 'MarkerEdgeColor', rightColor); % /10 to convert it to cm
                        % Plot a red line from the step point to the projection on the infinite line
                        plot([step_r(1, 1)/10, proj_point_r(1)/10], [step_r(1, 2)/10, proj_point_r(2)/10], 'r-', 'LineWidth', 2);
    
                        frame = getframe(gcf); % Get the current figure frame
                        writeVideo(videoWriter, frame); % Write the frame to the video
                 end
            
                % % Capture the current frame
                % frame = getframe(gcf); % Get the current figure frame
                % writeVideo(videoWriter, frame); % Write the frame to the video
                % 
                % % Pause to control the frame rate (1 second)
                % pause(pause_t); % Pause for 1 second before the next iteration
                  
                % Increment iteration counter
                iterationBuffer = iterationBuffer + 1; 
                if iterationBuffer == max_iteration
                    iterationBuffer = 1;
                    clf
                end
            end %BREAK HERE
            
            % Close the video file
            close(videoWriter);
        end
    end % end plot_longitudinal_error == 1
end



