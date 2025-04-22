Brief summary of the main scripts and the functions they use of this project to analyze motion capture data. If the functions are not listed here, they are not used for MoCap processing (they might be for EMG or IMU processing).:

RunMultipleTrials: this script must be used first to process the raw data from the mat files of each trial of each individual participants, and it generates an excel file with some parameters for visual inspection and another mat file with the calculated variables, which will be used by the analysis script.

Locatemarkers: locate all markers used
LocateSteps: detect the steps of each trial
StepError: calculate the absolute and signed error for the straights and curves. It is possible to generate a video of the step distribution and errors by setting the "video" variable to 1 and specifying the trial in line 10
StepParametersHuxham: calculate step length, step length var, step width, and step width var
WalkingSpeed: calculate walking speed and var, and stride time and var
Head_back_angles: calculate head and back angles
Dist_gaze: calculate how ahead participants are looking
LocateHeelStrikes: not used anymore
CalculateThreshold: not used anymore
Error_single_trial: not used anymore, but it served to analyse just a single trial from a participant

Multiple/Single _Participants: it plots the results in different ways for all participants together or only for one participant. Some features can be selected in lines 9.11, and also the type of plot to be made can be chosen. To know what each plot is like, read the comments in these scripts.

Sort_bylabel: one of the most important functions. It sorts the data in different ways that are used in the following functions
BaselineErrors: it calculates the headset baselines for each participant
BoxPlot_error
Accuracy_bars_error
Multiple_customized
Metabolics
All_conditions
Performance_Prompts
Compare_baselines
Survey_responses
Correlation_matrix
Correlation_prompts: used to calculate a few correlations
PerturbationsvsStepWidth: not used anymore

RunTrialData: not used anymore

SpeedPlotting: not used anymore