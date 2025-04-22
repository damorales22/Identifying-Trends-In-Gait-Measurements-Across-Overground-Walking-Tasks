%% Calculate head and trunk angles
function [head_angles, stj_headset_angles, trunk_angles] = Head_back_angles(trial, markerdata, trialname, static, bmh)
   
    %% Head angle (in the y-z plane)
    
    try
        coords_H1 = permute(markerdata.H1(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_H2 = permute(markerdata.H2(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_H3 = permute(markerdata.H3(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
    catch
        coords_H1 = permute(markerdata.HS1(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_H2 = permute(markerdata.HS1(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_H3 = permute(markerdata.HS3(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
    end

    delta_z = (coords_H2(:, 3) - coords_H3(:, 3))/10; % z-component difference; in cm
    d = sqrt(sum((coords_H3(1,:) - coords_H2(1,:)).^2))/10; % 3D distance between the markers (~15 cm)

    angle_radians = asin(delta_z / d);
    angle_degrees = rad2deg(angle_radians);

    angle_degrees = angle_degrees(~any(isnan(angle_degrees), 2), :); % Filter out rows with NaN

    head_angles = angle_degrees;
    
    if bmh == "BMH20" || bmh == "BMH17" || bmh == "BMH01" % In these trial, markers H2 and H3 are exchanged
        head_angles = -head_angles;
    end
    
    % figure
    % plot((1:11000)/120,head_angles(1:11000))

    %% Alternative qualitative method - STJ-headset angle
    try
        coords_STJ = permute(markerdata.STJ(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        delta_z = (coords_H3(:, 3) - coords_STJ(:, 3))/10; % z-component difference; in cm
        d = sqrt(sum((coords_STJ(1,:) - coords_H3(1,:)).^2))/10; % Distance (escalar) between the markers 
    
        angle_radians = asin(delta_z / d);
        angle_degrees = rad2deg(angle_radians);
    
        angle_degrees = angle_degrees(~any(isnan(angle_degrees), 2), :); % Filter out rows with NaN
        stj_headset_angles = real(angle_degrees);
    catch % for example, if there is no STJ marker
        stj_headset_angles = NaN;
    end
    
    % figure
    % plot((1:12000)/120,stj_headset_angles(1:12000))

    %% Trunk angle

    % To account for the movements or mismatch in the pelvis sensors, I do the average of left and right
    try
        coords_RPSIS = permute(markerdata.RPSIS(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_LPSIS = permute(markerdata.LPSIS(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z
        coords_C7 = permute(markerdata.C7(1,1:3,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y,z

        coords_pelvis = (coords_RPSIS + coords_LPSIS)/2;

        delta_z = (coords_pelvis(:, 3) - coords_C7(:, 3))/10; % z-component difference
        r = sqrt(sum((coords_C7 - coords_pelvis).^2, 2))/10;
    
        angle_radians = acos(delta_z ./ r);
        angle_degrees = rad2deg(angle_radians);
    
        angle_degrees = angle_degrees(~any(isnan(angle_degrees), 2), :); % Filter out rows with NaN
    
        trunk_angles = 180 - mean(angle_degrees);
    catch % if any of the markers is missing
        trunk_angles = NaN;
    end
    
    %figure
    %plot((1:14632)/120,trunk_angles)
   
end