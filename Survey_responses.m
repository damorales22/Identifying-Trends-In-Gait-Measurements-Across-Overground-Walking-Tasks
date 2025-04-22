%function [] = Survey_responses (imported_data,error,baseline_errors,metabolic_data,bmh)
function [] = Survey_responses (data_labels, plot_option, participants, show_cluster_analysis, show_formula, normalize_metabolics, bmh)

    if participants == "multiple"
        bmh_text = "Across all Participants";
    else
        bmh_text = string(bmh);
    end

    list_labels = fieldnames(data_labels); %labels stored in the struct

    %% Difficulty levels
    difficulty_levels = {'Extremely Easy', 'Somewhat Easy', 'Neither Easy nor Difficult', 'Somewhat Difficult', 'Extremely Difficult'};

    %% Colors and shapes
    balance_codes = zeros(1, 27);
    accuracy_codes = zeros(1, 27);
    speed_codes = zeros(1, 27);

    for i = 1:length(list_labels)
        trial = list_labels{i};
        
        % Extract the speed, accuracy, balance from the trial string (sxaybz)
        s = str2double(trial(2));  % Speed
        a = str2double(trial(4));  % Accuracy
        b = str2double(trial(6));  % Balance
        
        % According to the order stored in the struct
        accuracy_codes(i) = a; % Map accuracy to color (0=green, 1=red, 2=blue)
        balance_codes(i) = b; % Map balance to shape (0=Circle, 1=Square, 2=Triangle)  
        speed_codes(i) = s;
    end

    if plot_option == 1 % Color: accuracy. Shape: balance
        color_codes = accuracy_codes;
        shape_codes = balance_codes;
        shapes = {'o', 's', '^'};  % Circle, Square, Triangle (for balance)
        colors = {'g', 'r', 'b'};  % Green, Red, Blue (for accuracy)
        legend_labels = {'Accuracy: None, Disturbance: None', 'Accuracy: None, Disturbance: Medium', ...
                        'Accuracy: None, Disturbance: High', 'Accuracy: Medium, Disturbance: None', ...
                        'Accuracy: Medium, Disturbance: Medium', 'Accuracy: Medium, Disturbance: High', ...
                        'Accuracy: High, Disturbance: None', 'Accuracy: High, Disturbance: Medium', ...
                        'Accuracy: High, Disturbance: High'};
        clusters_color_silhouette = {'No Accuracy', 'Medium Accuracy', 'High Accuracy'};
    elseif plot_option == 2 % Color: balance. Shape: speed
        color_codes = balance_codes;
        shape_codes = speed_codes;
        shapes = {'p', 'd', '<'};  
        colors = {'#00B5B8', '#FF00FF', '#FFA500'};  % Cyan (#00B5B8), Magenta (#FF00FF), Orange (#FFA500)
        legend_labels = {'Disturbance: None, Speed: Slow', 'Disturbance: None, Speed: Medium', ...
                        'Disturbance: None, Speed: Fast', 'Disturbance: Medium, Speed: Slow', ...
                        'Disturbance: Medium, Speed: Medium', 'Disturbance: Medium, Speed: Fast', ...
                        'Disturbance: High, Speed: Slow', 'Disturbance: High, Speed: Medium', ...
                        'Disturbance: High, Speed: Fast'};
        clusters_color_silhouette = {'No Disturbance', 'Medium Disturbance', 'High Disturbance'};
    elseif plot_option == 3 % Color: speed. Shape: accuracy
        color_codes = speed_codes;
        shape_codes = accuracy_codes;
        shapes = {'h', 'o', 'v'};  
        colors = {'#DAA520', '#FF1493', '#32CD32'};  % Gold (#DAA520), DeepPink (#FF1493), LimeGreen (#32CD32)
        legend_labels = {'Speed: Slow, Accuracy: None', 'Speed: Slow, Accuracy: Medium', ...
                        'Speed: Slow, Accuracy: High', 'Speed: Medium, Accuracy: None', ...
                        'Speed: Medium, Accuracy: Medium', 'Speed: Medium, Accuracy: High', ...
                        'Speed: Fast, Accuracy: None', 'Speed: Fast, Accuracy: Medium', ...
                        'Speed: Fast, Accuracy: High'};
        clusters_color_silhouette = {'Slow Speed', 'Medium Speed', 'Fast Speed'};
    else
        error("Not valid plot option input")
    end

    %% Walking speed by speed survey (normalized) - Plot 1 
    % Prompt speed as color

    figure
    hold on
    
    x_all = [];
    y_all = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;
    
    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).walking_speed);
        y = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape 
        color = colors{color_codes(i) + 1};  % Get color 

        if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end

        %plot(x, y, 'o', 'Marker', shape, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 13, 'HandleVisibility', 'off');
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.01, max(x_all)+0.01, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all) - 0.07, max(y_all), formula_text, 'FontSize', 16);
    end

    %%legend show;
    title("Walking Speed by Perceived Walking Speed Priority (Normalized)");
    set(gca, 'FontSize', 20);
    xlabel('Walking Speed (m/s)');
    ylabel('Perceived Walking Speed Priority (Normalized)');
    xlim([min(x_all) - 0.02; max(x_all) + 0.02])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'best');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        % Silhouette - By accuracy across all conditions (color)
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end

    %% Step width variability by balance survey (normalized) - Plot 2
    % Applied balance disturbance as color

    figure
    hold on
    
    x_all = [];
    y_all = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};

        x = mean(data_labels.(current_label).step_width_var);
        %y = mean(data_labels.(current_label).surveys.balance_normalized);
        y = mean(data_labels.(current_label).survey_balance);

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        x_all = [x_all; x]; % Accumulate data for regression
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
 
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.1, max(x_all)+0.1, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        formula_text = sprintf('y = %.3fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        text(max(x_all) - 1, min(y_all) , formula_text, 'FontSize', 16);
    end
    
    %legend show;
    title("Step Width Variability by Perceived Balance Priority (Normalized)");
    set(gca, 'FontSize', 20);
    xlabel('Step Width Variability (cm)');
    ylabel('Perceived Balance Priority (Normalized)');
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    
    legend_handle.FontSize = 14; % size of the legend
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        % Silhouette - By accuracy across all conditions (color)
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end

    %% Actual metabolics by metabolics survey (normalized) - Plot 3

    %for i=1:3
    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;
    %points_a0b0 = []; points_a0b1 = []; points_a0b2 = []; points_a1b0 = []; points_a1b1 = []; points_a1b2 = []; points_a2b0 = []; points_a2b1 = []; points_a2b2 = []; 

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
        
        if normalize_metabolics == 0
            x = mean(data_labels.(current_label).metabolics.not_normalized);
        elseif normalize_metabolics == 1
            x = mean(data_labels.(current_label).metabolics.normalized);
        end
        
        y = mean(data_labels.(current_label).survey_metabolics);

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end

        x_all = [x_all; x]; % Accumulate data for regression
        y_all = [y_all; y];
        
        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
  
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   

    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all), max(x_all), 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        formula_text = sprintf('y = %.4fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        text(min(x_all)-0.3, min(y_all)+0.001, formula_text, 'FontSize', 16);
    end
    %text(475, 0.2, formula_text, 'FontSize', 12);

    set(gca, 'FontSize', 20);
    if normalize_metabolics == 0
        xlabel('Energy Expenditure (W)');
        xlim([min(x_all) - 5; max(x_all) + 5])
        ylim([min(y_all) - 0.01; max(y_all) + 0.01])
    elseif normalize_metabolics == 1
        xlabel('Energy Expenditure (W/kg)');
        xlim([min(x_all) - 0.4; max(x_all) + 0.4])
        ylim([min(y_all) - 0.005; max(y_all) + 0.005])
    end
    
    ylabel('Perceived Energy Expenditure Reduction Priority (Normalized)');
    title("Energy Expenditure by Perceived Energy Expenditure Reduction Priority");
    %legend show;
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries

    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end
    %end

    %% Accuracy (step error) by accuracy survey (normalized) - Plot 4
    % Prompt accuracy as color

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};

        x = mean(data_labels.(current_label).step_error);
        %y = mean(data_labels.(current_label).surveys.foot_placement_normalized);
        y = mean(data_labels.(current_label).survey_foot_placement);

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        x_all = [x_all; x]; % Accumulate data for regression
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   

    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.1, max(x_all)+0.1, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        formula_text = sprintf('y = %.4fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        text(max(x_all)-1.2, max(y_all)-0.01, formula_text, 'FontSize', 16);
    end
    
    %legend show;
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([min(y_all) - 0.01; max(y_all) + 0.01])
    title("Step Error by Perceived Foot Placement Priority (Normalized)");
    set(gca, 'FontSize', 20);
    xlabel('Step Error (cm)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'best');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end

    %% Step width var by difficulty survey (normalized) - Plot 5
    % Prompt speed as color

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_width_var);
        %y = mean(data_labels.(current_label).surveys.walking_speed_normalized);
        y = mean(data_labels.(current_label).survey_difficulty);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.1, max(x_all)+0.1, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-0.9, min(y_all)-0.6, formula_text, 'FontSize', 16);
    end

    % Show difficulty levels
    yticks(1:5);
    yticklabels(difficulty_levels);

    %legend show;
    %legend('Location', 'southeast');
    title("Step Width Variability by Perceived Difficulty");
    set(gca, 'FontSize', 20);
    xlabel('Step Width Variability (cm)');
    %ylabel('Perceived Difficulty');
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([1 5.05])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end

    %% Step error by difficulty survey (normalized) - Plot 6
    % Prompt speed as color

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_error);
        %y = mean(data_labels.(current_label).surveys.walking_speed_normalized);
        y = mean(data_labels.(current_label).survey_difficulty);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.3, max(x_all)+0.3, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all) -0.7 , max(y_all), formula_text, 'FontSize', 16);
    end

    % Show difficulty levels
    yticks(1:5);
    yticklabels(difficulty_levels);

    %legend show;
    % legend('Location', 'southeast');
    title("Step Error by Perceived Difficulty");
    set(gca, 'FontSize', 20);
    xlabel('Step Error (cm)');
    %ylabel('Perceived Difficulty');
    xlim([min(x_all) - 0.5; max(x_all) + 0.5])
    ylim([1 5.05])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'best');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_cluster_analysis == 1
        % Silhouette - By accuracy across all conditions (color)
        clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Accuracy (Color)');
        yticklabels(clusters_color_silhouette);

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Silhouette - By balance across all conditions (shape)
        clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
        labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);

        % Silhouette - By speed for that condition (shape + color)
        clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
        labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot by Speed for that Condition (Color + Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by accuracy (color)
        clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Accuracy (Color)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;

        % Mean by shape (balance)
        clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
        labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
        
        figure;
        silhouette(clusters, labels_silhouette);
        
        xlabel('Silhouette Value');
        ylabel('Cluster');
        title('Silhouette Plot of the Means by Balance (Shape)');

        silhouetteValues = silhouette(clusters, labels_silhouette);
        avgSilhouette = mean(silhouetteValues);
        disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
        ax = gca;               
        ax.FontSize = 20;
    end

    %% Walking speed priority by balance priority - Plot 7

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).survey_walking_speed);
        y = mean(data_labels.(current_label).survey_balance);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all), formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Perceived Balance Priority");
    set(gca, 'FontSize', 20);
    xlabel('Perceived Walking Speed Priority (Normalized)');
    ylabel('Perceived Balance Priority (Normalized)');
    xlim([min(x_all) - 0.01; max(x_all) + 0.01])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

    %% Balance prior by difficulty - Plot 8

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;
    
    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        y = mean(data_labels.(current_label).survey_difficulty);
        x = mean(data_labels.(current_label).survey_balance);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all) - 0.5, min(y_all)+0.005, formula_text, 'FontSize', 16);
    end

    yticks(1:5);
    yticklabels(difficulty_levels);

    %legend show;
    title("Perceived Difficulty by Perceived Balance Priority");
    set(gca, 'FontSize', 20);
    %xlabel('Perceived Difficulty');
    xlabel('Perceived Balance Priority (Normalized)');
    ylim([1; 5.05])
    xlim([min(x_all) - 0.01; max(x_all) + 0.01])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'southeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Walking speed prior by foot placement prior - Plot 9

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).survey_walking_speed);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.008, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Perceived Walking Speed Priority");
    set(gca, 'FontSize', 20);
    xlabel('Perceived Walking Speed Priority (Normalized)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 0.01; max(x_all) + 0.01])
    ylim([min(y_all) - 0.01; max(y_all) + 0.01])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Walking speed prior by difficulty - Plot 10

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        y = mean(data_labels.(current_label).survey_difficulty);
        x = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)-0.25, formula_text, 'FontSize', 16);
    end

    yticks(1:5);
    yticklabels(difficulty_levels);

    %legend show;
    title("Perceived Walking Speed Priority by Perceived Difficulty");
    set(gca, 'FontSize', 20);
    xlabel('Perceived Walking Speed Priority (Normalized)');
    %xlabel('Perceived Difficulty');
    ylim([1; 5.05])
    xlim([min(x_all) - 0.01; max(x_all) + 0.01])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% Foot placement priority by head angle - Plot 11 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).head_angle);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-6, min(y_all)-0.005, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Head Angle");
    set(gca, 'FontSize', 20);
    xlabel('Head Angle (deg)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 1; max(x_all) + 1])
    ylim([min(y_all) - 0.03; max(y_all) + 0.03])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Perceived Difficulty by head angle - Plot 12 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).head_angle);
        y = mean(data_labels.(current_label).survey_difficulty);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.004, formula_text, 'FontSize', 16);
    end

    yticks(1:5);
    yticklabels(difficulty_levels);
    title("Perceived Difficulty by Head Angle");
    set(gca, 'FontSize', 20);
    xlabel('Head Angle (deg)');
    %ylabel('Perceived Difficulty');
    ylim([1; 5.05])
    xlim([min(x_all) - 2; max(x_all) + 2])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% Perceived Difficulty by trunk angle - Plot 13

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).trunk_angle);
        y = mean(data_labels.(current_label).survey_difficulty);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.004, formula_text, 'FontSize', 16);
    end

    yticks(1:5);
    yticklabels(difficulty_levels);
    title("Perceived Difficulty by Trunk Angle");
    set(gca, 'FontSize', 20);
    xlabel('Trunk Angle (deg)');
    %ylabel('Perceived Difficulty');
    ylim([1; 5.05])
    xlim([min(x_all) - 0.4; max(x_all) + 0.4])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% Walking speed priority by step length var - Plot 14

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_length_var);
        y = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all)-0.5, min(y_all)-0.12, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Step Length Variability");
    set(gca, 'FontSize', 20);
    xlabel('Step Length Variability (cm)');
    ylabel('Perceived Walking Speed Priority (Normalized)');
    xlim([min(x_all) - 1; max(x_all) + 1])
    ylim([min(y_all) - 0.15; max(y_all) + 0.05])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% Foot placement priority by step error var - Plot 15 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_error_std);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.004, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Step Error Variability");
    set(gca, 'FontSize', 20);
    xlabel('Step Error Variability (cm)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;


%% Walking speed priority by head angle - Plot 16

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).head_angle);
        y = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.004, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Head Angle");
    set(gca, 'FontSize', 20);
    xlabel('Head Angle (deg)');
    ylabel('Perceived Walking Speed Priority (Normalized)');
    xlim([min(x_all) - 1.5; max(x_all) + 1.5])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

    %% Foot placement priority by head angle - Plot 17

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).head_angle);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-5, min(y_all), formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Head Angle");
    set(gca, 'FontSize', 20);
    xlabel('Head Angle (deg)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 1.5; max(x_all) + 1.5])
    ylim([min(y_all) - 0.03; max(y_all) + 0.03])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% EE reduction priority by walking speed prior - Plot 18

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).survey_walking_speed);
        y = mean(data_labels.(current_label).survey_metabolics);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-0.05, min(y_all), formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Perceived Energy Expenditure Reduction Prioirity");
    set(gca, 'FontSize', 20);
    ylabel('Perceived Energy Expenditure Reduction Prioirity (Normalized)');
    xlabel('Perceived Walking Speed Priority (Normalized)');
    %xlim([min(x_all) - 1.5; max(x_all) + 1.5])
    ylim([0.06 0.17])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

    %% Foot placement priority by step width - Plot 19 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_width);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all), 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-1.2, min(y_all)-0.005, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Step Width");
    set(gca, 'FontSize', 20);
    xlabel('Step Width (cm)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([min(y_all) - 0.03; max(y_all) + 0.03])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Foot placement priority by step error var - Plot 20 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_error);
        y = mean(data_labels.(current_label).survey_foot_placement);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)+0.004, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Step Error");
    set(gca, 'FontSize', 20);
    xlabel('Step Error (cm)');
    ylabel('Perceived Foot Placement Priority (Normalized)');
    xlim([min(x_all) - 0.4; max(x_all) + 0.4])
    ylim([min(y_all) - 0.02; max(y_all) + 0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% Walking speed priority by step width - Plot 21 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).step_width);
        y = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all), 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all)-0.005, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Step Width");
    set(gca, 'FontSize', 20);
    xlabel('Step Width (cm)');
    ylabel('Perceived Walking Speed Priority (Normalized)');
    xlim([min(x_all) - 0.2; max(x_all) + 0.2])
    ylim([min(y_all) - 0.03; max(y_all) + 0.03])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% EE priority by head angle - Plot 22 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).head_angle);
        y = mean(data_labels.(current_label).survey_metabolics);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all), min(y_all), formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Energy Expenditure Reduction Priority by Head Angle");
    set(gca, 'FontSize', 20);
    xlabel('Head Angle (deg)');
    ylabel('Energy Expenditure Reduction Priority (Normalized)');
    xlim([min(x_all) - 1; max(x_all) + 1])
    ylim([min(y_all) - 0.005; max(y_all) + 0.005])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Walking speed priority by speed std - Plot 23 

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).walking_speed_std, 'omitnan');
        y = mean(data_labels.(current_label).survey_walking_speed);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all), max(x_all), 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all)-0.003, min(y_all)+0.01, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Walking Speed Priority by Walking Speed Variability");
    set(gca, 'FontSize', 20);
    xlabel('Walking Speed Variability (m/s)');
    ylabel('Walking Speed Priority (Normalized)');
    xlim([min(x_all)-0.006; max(x_all)+0.008])
    ylim([min(y_all)-0.02; max(y_all)+0.02])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    
    legend_handle.FontSize = 14; % size of the legend

    hold off;

%% Metabolics by difficulty - Plot 24

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        y = mean(data_labels.(current_label).survey_difficulty);
        x = mean(data_labels.(current_label).metabolics.normalized);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(max(x_all)-0.5, min(y_all)-0.25, formula_text, 'FontSize', 16);
    end

    yticks(1:5);
    yticklabels(difficulty_levels);

    %legend show;
    title("Energy Expenditure by Perceived Difficulty");
    set(gca, 'FontSize', 20);
    xlabel('Energy Expenditure (W/kg)');
    %xlabel('Perceived Difficulty');
    ylim([1; 5.05])
    xlim([min(x_all) - 0.1; max(x_all) + 0.1])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

%% EE reduction priority by foot placement prior - Plot 25

    figure
    hold on
    
    x_all = [];
    y_all = [];

    points.low = [];
    points.medium = [];
    points.high = [];

    count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1; count7 = 1; count8 = 1; count9 = 1;

    for i = 1:size(list_labels, 1)
        current_label = list_labels{i};
                
        x = mean(data_labels.(current_label).survey_foot_placement);
        y = mean(data_labels.(current_label).survey_metabolics);
        x_all = [x_all; x]; % Accumulate data for regression line
        y_all = [y_all; y];

        shape = shapes{shape_codes(i) + 1};  % Get shape for balance
        color = colors{color_codes(i) + 1};  % Get color for accuracy

         if shape_codes(i) == 0 && color_codes(i) == 0
            points_a0b0(count1,1) = x;
            points_a0b0(count1,2) = y;
            count1 = count1+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 0
            points_a0b1(count2,1) = x;
            points_a0b1(count2,2) = y;
            count2 = count2+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 0
            points_a0b2(count3,1) = x;
            points_a0b2(count3,2) = y;
            count3 = count3+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 1
            points_a1b0(count4,1) = x;
            points_a1b0(count4,2) = y;
            count4 = count4+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 1
            points_a1b1(count5,1) = x;
            points_a1b1(count5,2) = y;
            count5 = count5+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 1
            points_a1b2(count6,1) = x;
            points_a1b2(count6,2) = y;
            count6 = count6+1;
        elseif shape_codes(i) == 0 && color_codes(i) == 2
            points_a2b0(count7,1) = x;
            points_a2b0(count7,2) = y;
            count7 = count7+1;
        elseif shape_codes(i) == 1 && color_codes(i) == 2
            points_a2b1(count8,1) = x;
            points_a2b1(count8,2) = y;
            count8 = count8+1;
        elseif shape_codes(i) == 2 && color_codes(i) == 2
            points_a2b2(count9,1) = x;
            points_a2b2(count9,2) = y;
            count9 = count9+1;
        end
        
        scatter(x, y, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.1);   % Reduce opacity
    end
    
    % Plot averages
    avg_a0b0 = mean(points_a0b0,1);
    scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b1 = mean(points_a0b1,1);
    scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a0b2 = mean(points_a0b2,1);
    scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});

    avg_a1b0 = mean(points_a1b0,1);
    scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b1 = mean(points_a1b1,1);
    scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a1b2 = mean(points_a1b2,1);
    scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});

    avg_a2b0 = mean(points_a2b0,1);
    scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b1 = mean(points_a2b1,1);
    scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});

    avg_a2b2 = mean(points_a2b2,1);
    scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});   
    
    % Perform linear regression
    coefficients = polyfit(x_all, y_all, 1); % Fit line: y = mx + c
    y_fit = polyval(coefficients, x_all);   % Predicted y values
    m = coefficients(1); % Slope
    c = coefficients(2); % Intercept
    
    % Compute R^2
    ss_total = sum((y_all - mean(y_all)).^2); % Total sum of squares
    ss_residual = sum((y_all - y_fit).^2);   % Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total);
    
    % Plot the regression line
    x_range = linspace(min(x_all)-0.005, max(x_all)+0.005, 100); % Create range for the line
    y_range = polyval(coefficients, x_range);
    plot(x_range, y_range, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    
    % Display the formula and R^2 in a location away from the line
    if show_formula == 1
        if c >= 0
            formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
        else
            formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
        end
        text(min(x_all)-0.03, min(y_all)-0.002, formula_text, 'FontSize', 16);
    end

    %legend show;
    title("Perceived Foot Placement Priority by Perceived Energy Expenditure Reduction Prioirity");
    set(gca, 'FontSize', 20);
    ylabel('Perceived Energy Expenditure Reduction Prioirity (Normalized)');
    xlabel('Perceived Foot Placement Priority (Normalized)');
    %xlim([min(x_all) - 1.5; max(x_all) + 1.5])
    ylim([0.06 0.17])
    grid on;

    legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
    
    index = 1;
    for color = colors
        for shape = shapes
            legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
            index = index + 1;
        end
    end

    legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
    legend_handle.FontSize = 14; % size of the legend
    hold off;

end