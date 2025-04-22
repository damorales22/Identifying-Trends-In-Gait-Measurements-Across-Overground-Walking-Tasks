# import modules
import csv
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


## load heel strike frames

# define file names
strikes_r_path = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/strikes_r.csv'
strikes_l_path = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/strikes_l.csv'

# .csv loading function
def loadstrikes(filename):
    data_dict = {}

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        
        # Read all rows from the CSV file
        rows = list(csvreader)
        
        # Assuming there are exactly 9 columns in each row
        num_trials = len(rows[0])  # number of trials (columns)
        
        for i in range(num_trials):
            trial_name = f"Trial{i+1:04}"  # Generate trial names like Trial0001, Trial0002, etc.
            data_values = [float(row[i]) for row in rows]  # Collect data from the i-th column
            data_values = [value for value in data_values if not math.isnan(value)] # Removes nans from end of trial data
            data_dict[trial_name] = data_values
    
    return data_dict

# load strikes into dictionaries
strikes_r = loadstrikes(strikes_r_path)
strikes_l = loadstrikes(strikes_l_path)


## load IK data

# .csv loading function
def loadtrials(directory):
    trial_dicts = {}

    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            trial_name = os.path.splitext(filename)[0]  # Extract trial name from filename
            trial_dict = {}

            # Open the CSV file
            with open(os.path.join(directory, filename), 'r') as csvfile:
                csvreader = csv.reader(csvfile)
                headers = next(csvreader)  # Read the headers from the first row
                
                # Initialize lists for each header
                columns = {header: [] for header in headers}

                # Read the data rows
                for row in csvreader:
                    for idx, value in enumerate(row):
                        columns[headers[idx]].append(float(value))  # Assuming all values are numeric

            # Store the columns in the trial_dict
            for header, values in columns.items():
                trial_dict[header] = values
            
            # Store the trial_dict in the trial_dicts with the trial_name as key
            trial_dicts[trial_name] = trial_dict

    return trial_dicts

# load IK data into trials dictionary
directory = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/'  # Replace with the directory containing your CSV files
trials = loadtrials(directory)

## Segmentation and Plotting Function
def AnglePlot(trialname, strikes, label, subplotname):
    # Get angle data
    data = np.array(trials[trialname][label])

    # Convert to degrees
    data = list(data * (180/np.pi))

    ## Segment Angle Data
    num_segments = len(strikes) - 1 # don't include last partial step
    segments = []

    # Segment the data based on time stamps
    for i in range(num_segments):
        start_ind = strikes[i] # start index of segment
        end_ind = strikes[i+1] # end index of segment
        segments.append(data[int(start_ind-1):int(end_ind)]) # Appends segment onto list

    # Create steps (based on longest segment)
    max_length = len(max(segments,key=len)) # find longest segment
    steps = np.linspace(0,100, max_length, endpoint=False)

    # Initialize arrays for average y data
    avg = np.zeros(len(steps))

    # Initialize a counter for valid segments
    valid_segments_count = 0
    # Calculate average of all segment amplitudes
    total_amp = 0 # Accumulate all segment amplitudes
    for i in range(num_segments):
        tempmax = max(segments[i])
        total_amp += tempmax
    average_amplitude = total_amp/num_segments # avg amplitude of all segments

    # Define a threshold based on avg amplitude
    threshold2 = average_amplitude*3

    # Plot the segments and calculate average y data
    for i in range (num_segments):
        # Create stepping vector
        segment_length = len(segments[i])
        percent = np.linspace(0,100, segment_length, endpoint=False)

        # Check if segment should be included
        if any(segments[i] > threshold2):
            continue
        # Plot the segment
        subplotname.plot(percent, segments[i], color='#c0c0c0')

        # Interpolate amplitude to match steps
        interpolator = interp1d(percent, segments[i], kind='linear', fill_value='extrapolate')
        interpolated = np.array(interpolator(steps))

        # Accumulate average y data
        avg = np.add(avg, interpolated)

        # Increment count of valid segements
        valid_segments_count += 1

    # Calculate average by dividing
    if valid_segments_count > 0:
        avg = avg / valid_segments_count
    
    # Plot the avg y data
    subplotname.plot(steps, avg, color='#003366', linewidth=2)

    # Plot Appearance
    subplotname.set_xlabel('Gait Cycle (%)')
    subplotname.set_ylabel('Joint Angle (deg)')
    subplotname.axis('tight')
    subplotname.set_title("Trial " + trialname[len(trialname)-1])

# Variables
trialtitle = "Sarah's Trials - "
trialprefix = "Trial000"
numtrials = 9
angles = ["ankle_angle", "hip_flexion", "knee_angle"]
joints = ["Ankle", "Hip", "Knee"]

# Creating Angle Labels
labels = []
jointlabels = []
for i in range(len(angles)):
    labels.append(angles[i] + '_r')
    jointlabels.append(' ' + joints[i] + ' Angle')
for i in range(len(angles)):
    labels.append(angles[i] + '_l')
    jointlabels.append(' ' + joints[i] + ' Angle')

# Loop through labels
figures = {}
for i in range(len(labels)):
    label = labels[i]
    jointlabel = jointlabels[i]
    
    # Check foot
    if label.endswith("r"):
        side = "right"
    else:
        side = "left"
    
    # Create figure
    figures[str(i)] = plt.figure(figsize=(12,7))
    subplotsdict = {}
    
    # Gets data per trial
    for j in range(1, numtrials+1):
        trialname = trialprefix + str(j)
        
        # Save strides from desired foot
        if side == "right":
            strikes = strikes_r[trialname]
        else:
            strikes = strikes_l[trialname]
        
        # Add new subplot to the graph
        subplotsdict[trialname] = figures[str(i)].add_subplot(3, int(numtrials/3), j)
        
        # Run Plotting function
        AnglePlot(trialname, strikes, label, subplotsdict[trialname])

    # Sets Figure Title    
    figures[str(i)].suptitle(trialtitle + side.replace('r', 'R').replace('l', 'L') + jointlabel)

    # Loop through all axes to mins and maxs
    allYLim = np.zeros([numtrials,2])
    for j in range(numtrials):
        trialname = trialprefix + str(j+1)
        allYLim[j, :] = subplotsdict[trialname].get_ylim() # Get current y-limits for each subplot (min, max)
    
    # Find the overall min and max
    maxYLim = max(allYLim[:,1])
    minYLim = min(allYLim[:,0])

    # Set standard y-limits
    for j in range(numtrials):
        trialname = trialprefix + str(j+1)
        subplotsdict[trialname].set_ylim(minYLim, maxYLim)
    
    # Spaces out graph
    figures[str(i)].tight_layout()

# Show all figures
plt.show()