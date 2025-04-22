%% Batchc3dExport.m
% Runs c3dExport.m for a set of Trials, saving .mot and .trc files for each
close all
clear
clc

%% Variables
filepath = "C:\Ability Lab\7_31_2024\";
input1 = 1; % use center of pressure
input2 = 0; % don't convert units
numtrials = 27;

%% Run c3dExport.m for all trials 
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    export(input1, input2, filepath, trialname)
end

function export(useCenterOfPressureAsMomentsPoint, convertLengthUnits, path, filename)
    % The function names output .trc and .mot files with the same
% basename as the selcted input .c3d file.

    
    if nargin < 2
        convertLengthUnits =  0
    end
    %% Example of using the Matlab-OpenSim class 
    
    %% Load OpenSim libs
    import org.opensim.modeling.*
    
    %% Get the path to a C3D file
    c3dpath = fullfile(path + filename  + '.c3d');
    disp('loaded' + c3dpath)
    
    %% Construct an opensimC3D object with input c3d path
    % Constructor takes full path to c3d file and an integer for forceplate
    % representation in output forces file (0 = electrical center, 1 = COP). 
    c3d = osimC3D(c3dpath,useCenterOfPressureAsMomentsPoint);
    
    %% Get some stats...
    % Get the number of marker trajectories
    nTrajectories = c3d.getNumTrajectories();
    % Get the marker data rate
    rMakers = c3d.getRate_marker();
    % Get the number of forces 
    nForces = c3d.getNumForces();
    % Get the force data rate
    rForces = c3d.getRate_force();
    
    % Get Start and end time
    t0 = c3d.getStartTime();
    tn = c3d.getEndTime();
    
    %% Rotate the data 
    c3d.rotateData('x',-90)
    
    %% Get the c3d in different forms
    % Get OpenSim tables
    markerTable = c3d.getTable_markers();
    forceTable = c3d.getTable_forces();
    % Get as Matlab Structures
    [markerStruct forceStruct] = c3d.getAsStructs();
    
    %% Convert COP (mm to m) and Moments (Nmm to Nm)
    if convertLengthUnits
        c3d.convertMillimeters2Meters();
    end
    
    %% Write the marker and force data to file
    % Define output file names
    basename = strtok(filename,'.');
    markersFilename = strcat(basename,'_markers.trc');
    
    switch useCenterOfPressureAsMomentsPoint
        case 0
            forcesFilename = strcat(basename,'_forces_EC.mot');
        case 1
            forcesFilename = strcat(basename,'_forces_COP.mot');
    end
    
    % Write marker data to trc file.
    % c3d.writeTRC()                       Write to dir of input c3d.
    % c3d.writeTRC('Walking.trc')          Write to dir of input c3d with defined file name.
    % c3d.writeTRC('C:/data/Walking.trc')  Write to defined path input path.
    c3d.writeTRC(markersFilename);
    
    % Write force data to mot file.
    % c3d.writeMOT()                       Write to dir of input c3d.
    % c3d.writeMOT('Walking.mot')          Write to dir of input c3d with defined file name.
    % c3d.writeMOT('C:/data/Walking.mot')  Write to defined path input path.
    c3d.writeMOT(forcesFilename);

end