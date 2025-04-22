%% This Script Combines .mot files into a single file
clear
close all
clc

%% Variables
filepath = "C:/Ability Lab/6_24_2024/OpenSim Files/Jordan/IK/";
trialprefix = "Trial000";
trialsuffix1 = "_markers_segment_";
trialsuffix2 = "_ik";
numtrials = 9;

%% Loop through trials and segments
for i = 1:numtrials % loop through trials

    % get trial
    trialname = trialprefix + num2str(i);

    % load first segment
    segment = 0;
    segmentname = trialsuffix1 + num2str(segment) + trialsuffix2;
    file = filepath + trialname + segmentname + num2str('.mot');
    data = importdata(file);
    
    % get header text
    headers = data.textdata;

    % get numerical data
    data = data.data;

    % go to second segment
    segment = 1;
    segmentname = trialsuffix1 + num2str(segment) + trialsuffix2;
    file = filepath + trialname + segmentname + num2str('.mot');

    while isfile(file)
        tempfile = importdata(file); % load next segment
        tempfile = tempfile.data; % remove headers 
        data = cat(1, data, tempfile); % add to previous segment(s)

        % go to next segment
        segment = segment + 1; 
        segmentname = trialsuffix1 + num2str(segment) + trialsuffix2;
        file = filepath + trialname + segmentname + num2str('.mot');
    end
    
    % update row count in data
    rowcount = size(data, 1);

    for j = 1:size(headers, 1)
        for k = 1: size(headers, 2)
            if ischar(headers{j, k}) && contains(headers{j, k}, 'nRows=')
                headers{j, k} = sprintf('nRows = %d', rowcount);
                break;  % Exit the loop once replacement is done
            end
        end
    end

    % Check if the file exists
    filename = filepath + trialname + '.mot'; % create filename

    if exist(filename, 'file') == 2
        % File exists
        choice = questdlg('The file already exists. Do you want to delete it and continue?', ...
                          'File Exists', ...
                          'Delete and Continue', 'Cancel', 'Cancel');
        if strcmp(choice, 'Delete and Continue')
            % Delete the existing file
            delete(filename);
            fprintf('Existing file deleted.\n');
        else
            fprintf('Operation cancelled.\n');
            return;  % Exit the script or function
        end
    end
    
    % Open file for writing
    fid = fopen(filename', 'w');

    % Write each row to the file
    for j = 1:size(headers, 1)
        for k = 1:size(headers, 2)
            fprintf(fid, '%s\t', headers{j,k});  % Assuming tab-delimited format
        end
        fprintf(fid, '\n');  % End of row
        % Check if we need to add an empty row after the 6th and 9th rows
        if j == 5 || j == 7
            fprintf(fid, '\n');  % Empty row
        end
    end
    
    % Write numerical data to .mot file
    for j = 1:size(data, 1)
        fprintf(fid, '%g\t', data(j, :));
        fprintf(fid, '\n');
    end
    
    % Close file
    fclose(fid);

    disp(trialname + ' combined')
end
