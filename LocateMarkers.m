% Locate Markers Function 
% function assumes that the floor markers are named 'FL1a', 'FL1b', 'FL2a',
% etc. with a markers being on the inside of the track and b markers being
% on the outside.

function [points] = LocateMarkers(markerdata,floordata)
    %% Remove non Floor markers
    
    % get field names
    %fieldNames = fieldnames(markerdata);
    fieldNames = fieldnames(floordata);
    
    % Initialize a logical index to track fields to keep
    keepField = false(size(fieldNames));
    
    % Identify fields containing 'FL'
    for i = 1:numel(fieldNames)
        if startsWith(fieldNames{i}, 'FL')    %we only wanna detect the floor (FL) markers
            keepField(i) = true;              
        end
    end
    
    %markerdata = rmfield(markerdata, fieldNames(~keepField)); %remove the fields that are not floor markers
    floordata = rmfield(floordata, fieldNames(~keepField)); %remove the fields that are not floor markers
    
    %% find Coordinates
    
    % initialize variables
    points = struct();
    coords = zeros(3, 3);
    
    %for i = 1:length(fields(markerdata))/2
    for i = 1:length(fields(floordata))/2 %Number of floor markers; /2 because there are 2 per tape
        %coord_i = mean(squeeze(markerdata.("FL" + num2str(i) + "a")(1, 1:3, :)), 2); % get INNER marker loc. 1:3 to take only the 3 forst rows; what is the fourth one? mean (X,2): mean accross the columns. squeeze to remove the first dimension
        %coord_o = mean(squeeze(markerdata.("FL" + num2str(i) + "b")(1, 1:3, :)), 2); % get OUTER marker loc
        %coord_m = (coord_o(:) + coord_i(:)).'/2; % gets the midpoint of the two (middle of the red tape)
        
        coord_i = mean(squeeze(floordata.("FL" + num2str(i) + "a")(1, 1:3, :)), 2); %first column of points
        coord_o = mean(squeeze(floordata.("FL" + num2str(i) + "b")(1, 1:3, :)), 2); %second column of points
        coord_m = (coord_o(:) + coord_i(:)).'/2;                                    %third column of points
 
        % Assign coordinates to the matrix
        coords(:, 1) = coord_i; %coords inner (x,y,z)
        coords(:, 2) = coord_o; %coords outer (x,y,z)
        coords(:, 3) = coord_m; %coords middle (x,y,z)
        points.("FL" + num2str(i)) = coords(:, :); %this is the output. 
    end
end