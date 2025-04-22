function [] = Compare_baselines(imported_data,error,baseline_errors,bmh)

    if nargin < 5
        bmh_text = "Across all Participants";
        participants = "multiple";
    else
        participants = "single";
        bmh_text = string(bmh);
    end

    %% Gather the data that can be plotted

    walking_speed = [];
    step_width = [];
    step_length = [];
    step_error = [];

    if participants == "single"
        walking_speed = imported_data.(bmh).walkingspeed;
        step_width = imported_data.(bmh).meanwidthstraights;
        step_width_var = imported_data.(bmh).straightsvariability;  
        step_length = imported_data.(bmh).meanlength;
        step_error = imported_data.(bmh).overall_errors.straight_mean_absolute;
    elseif participants == "multiple"
        for i=1:size(fieldnames(imported_data),1)
            names = string(fieldnames(imported_data));
            bmh = names(i);
            walking_speed(i,:) = imported_data.(bmh).walkingspeed;
            step_width(i,:) = imported_data.(bmh).meanwidthstraights;
            step_width_var(i,:) = imported_data.(bmh).straightsvariability;
            step_length(i,:) = imported_data.(bmh).meanlength_straights;
            step_error(i,:) = imported_data.(bmh).overall_errors.straight_mean_absolute;
        end
   end
   
    %% Classify trials according to their labels

    data_labels.s0a0b0 = [];
    data_labels.s0a0b1 = [];
    data_labels.s0a0b2 = [];
    data_labels.s0a1b0 = [];
    data_labels.s0a1b1 = [];
    data_labels.s0a1b2 = [];
    data_labels.s0a2b0 = [];
    data_labels.s0a2b1 = [];
    data_labels.s0a2b2 = [];
    
    data_labels.s1a0b0 = [];
    data_labels.s1a0b1 = [];
    data_labels.s1a0b2 = [];
    data_labels.s1a1b0 = [];
    data_labels.s1a1b1 = [];
    data_labels.s1a1b2 = [];
    data_labels.s1a2b0 = [];
    data_labels.s1a2b1 = [];
    data_labels.s1a2b2 = [];
    
    data_labels.s2a0b0 = [];
    data_labels.s2a0b1 = [];
    data_labels.s2a0b2 = [];
    data_labels.s2a1b0 = [];
    data_labels.s2a1b1 = [];
    data_labels.s2a1b2 = [];
    data_labels.s2a2b0 = [];
    data_labels.s2a2b1 = [];
    data_labels.s2a2b2 = [];

    list_labels = fieldnames(data_labels); %labels stored in the struct

    bmh_list = string(fieldnames(baseline_errors));

    for i=1:size(fieldnames(imported_data),1) %through bmhs
        names = string(fieldnames(imported_data));
        bmh = names(i);
        numtrials = size(imported_data.(bmh).labels,2);
        labels = imported_data.(bmh).labels;  % Format: 'sxaybz'

        for j = 3:29 %through trials
            current_label = labels{j}; % from the ordered label array
            trialname = sprintf('Trial%04d', j); 
            for k=1:size(list_labels,1)
                label_struct = list_labels{k};
                if strcmp(current_label, label_struct)
                    data_labels.(current_label).walking_speed(i,:) = walking_speed(i,j);
                    data_labels.(current_label).step_width(i,:) = step_width(i,j);
                    data_labels.(current_label).step_width_var(i,:) = step_width_var(i,j);
                    data_labels.(current_label).step_length(i,:) = step_length(i,j);
                    data_labels.(current_label).step_error(i,:) = step_error(i,j);
                else
                    continue
                end
            end
        end
    end

    % Initialize arrays for shape and color codes
    shape_codes = zeros(1, numtrials);
    color_codes = zeros(1, numtrials);

    for i = 1:length(list_labels)
        trial = list_labels{i};
        
        % Extract the speed, accuracy, balance from the trial string (sxaybz)
        s = str2double(trial(2));  % Speed
        a = str2double(trial(4));  % Accuracy
        b = str2double(trial(6));  % Balance
        
        % According to the order stored in the struct
        color_codes(i) = a; % Map accuracy to color (0=green, 1=red, 2=blue)
        shape_codes(i) = b; % Map balance to shape (0=Circle, 1=Square, 2yy=Triangle)  
    end

    %% Actual comparison

    fileID = fopen('C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Baseline check.txt', 'w'); % 'w' mode overwrites the file
    written = 0;
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        labels = imported_data.(bmh).labels;  % Format: 'sxaybz'
    
        % Loop for the s0a2b0 baseline
        for j = 3:29
            current_label = labels{j}; % from the ordered label array
            if double(imported_data.(bmh).prompt.speed(j,2)) == 0 && (double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 || double(imported_data.(bmh).prompt.accuracy(j,2)) == 1)
                if baseline_errors.(bmh).speed_s0a2b0 > data_labels.(current_label).step_error(i)
                    sentence = ['Baseline ' bmh '_s0a2b0: ' num2str(baseline_errors.(bmh).speed_s0a2b0) ...
                                '. Error in ' bmh '_' string(current_label) ': ' ...
                                num2str(data_labels.(string(current_label)).step_error(i)) '.'];
    
                    % Write the sentence to the file
                    fprintf(fileID, '%s', sentence);
                    fprintf(fileID, '\n');
                    written = 1;
                end
            end
        end
        if written == 1
            fprintf(fileID, '\n');
            written = 0;
        end
    
        % Loop for the s1a2b0 baseline
        for j = 3:29
            current_label = labels{j}; % from the ordered label array
            if double(imported_data.(bmh).prompt.speed(j,2)) == 1 && (double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 || double(imported_data.(bmh).prompt.accuracy(j,2)) == 1)
                if baseline_errors.(bmh).speed_s1a2b0 > data_labels.(current_label).step_error(i)
                    sentence = ['Baseline ' bmh '_s1a2b0: ' num2str(baseline_errors.(bmh).speed_s1a2b0) ...
                                '. Error in ' bmh '_' string(current_label) ': ' ...
                                num2str(data_labels.(string(current_label)).step_error(i)) '.'];
    
                    % Write the sentence to the file
                    fprintf(fileID, '%s', sentence);
                    fprintf(fileID, '\n');
                    written = 1;
                end
            end
        end
        if written == 1
            fprintf(fileID, '\n');
            written = 0;
        end
    
        % Loop for the s2a2b0 baseline
        for j = 3:29
            current_label = labels{j}; % from the ordered label array
            if double(imported_data.(bmh).prompt.speed(j,2)) == 2 && (double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 || double(imported_data.(bmh).prompt.accuracy(j,2)) == 1)
                if baseline_errors.(bmh).speed_s2a2b0 > data_labels.(current_label).step_error(i)
                    sentence = ['Baseline ' bmh '_s2a2b0: ' num2str(baseline_errors.(bmh).speed_s2a2b0) ...
                                '. Error in ' bmh '_' string(current_label) ': ' ...
                                num2str(data_labels.(string(current_label)).step_error(i)) '.'];
    
                    % Write the sentence to the file
                    fprintf(fileID, '%s', sentence);
                    fprintf(fileID, '\n');
                    written = 1;
                end
            end
        end
        if written == 1
            fprintf(fileID, '\n');
            written = 0;
        end
        fprintf(fileID, '\n');
    end
    
    % Close the file
    fclose(fileID);

end