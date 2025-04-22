# import modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

## load heel strike frames
# define file names
strikes_r_path = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/strikes_r.csv'
strikes_l_path = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/strikes_l.csv'

# .csv loading function
def load_strikes(filename):
    df = pd.read_csv(filename, header=None)
    data_dict = {f"Trial{i+1:04}": df[i].dropna().values for i in range(df.shape[1])}
    return data_dict

# load strikes into dictionaries
strikes_r = load_strikes(strikes_r_path)
strikes_l = load_strikes(strikes_l_path)

## load IK data
# .csv loading function
def load_trials(directory):
    trial_dicts = {}
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            trial_name = os.path.splitext(filename)[0]
            trial_dict = pd.read_csv(os.path.join(directory, filename)).to_dict(orient='list')
            trial_dicts[trial_name] = {key: np.array(value, dtype=float) for key, value in trial_dict.items()}
    return trial_dicts

# load IK data into trials dictionary
directory = '/Users/sarahli/Documents/Ability Lab/6_25_Data/Opensim Kinematics/'  # Replace with the directory containing your CSV files
trials = load_trials(directory)

## Segmentation and Plotting Function
def angle_plot(trialname, strikes, label, ax):
    # Get angle data and convert to degrees
    data = np.degrees(np.array(trials[trialname][label]))

    # Segment Angle Data
    num_segments = len(strikes) - 1
    segments = [data[int(strikes[i]-1):int(strikes[i+1])] for i in range(num_segments)]

    # Create steps (based on longest segment)
    max_length = max(len(segment) for segment in segments)
    steps = np.linspace(0, 100, max_length, endpoint=False)

    # Calculate average amplitude for threshold
    avg_amplitude = np.mean([max(segment) for segment in segments])
    threshold2 = avg_amplitude * 3

    # Initialize arrays for average y data
    avg = np.zeros(len(steps))
    valid_segments_count = 0

    # Plot the segments and calculate average y data
    for segment in segments:
        segment_length = len(segment)
        percent = np.linspace(0, 100, segment_length, endpoint=False)

        if np.any(segment > threshold2):
            continue

        # Interpolate amplitude to match steps
        interpolator = interp1d(percent, segment, kind='linear', fill_value='extrapolate')
        interpolated = interpolator(steps)

        # Accumulate average y data
        avg += interpolated
        valid_segments_count += 1
        ax.plot(percent, segment, color='#c0c0c0')

    # Calculate and plot average
    if valid_segments_count > 0:
        avg /= valid_segments_count
    ax.plot(steps, avg, color='#003366', linewidth=2)

    # Plot appearance
    ax.set_xlabel('Gait Cycle (%)')
    ax.set_ylabel('Joint Angle (deg)')
    ax.axis('tight')
    ax.set_title(f"Trial {trialname[-1]}")

# Variables
trialtitle = "Sarah's Trials - "
trialprefix = "Trial000"
numtrials = 9
angles = ["ankle_angle", "hip_flexion", "knee_angle"]
joints = ["Ankle", "Hip", "Knee"]

# Creating Angle Labels
labels = [f"{angle}_r" for angle in angles] + [f"{angle}_l" for angle in angles]
jointlabels = [' ' + joint + ' Angle' for joint in joints] * 2

# Loop through labels
figures = {}
for i, (label, jointlabel) in enumerate(zip(labels, jointlabels)):
    side = "right" if label.endswith("r") else "left"

    figures[str(i)] = plt.figure(figsize=(12, 7))
    subplotsdict = {}

    for j in range(1, numtrials + 1):
        trialname = f"{trialprefix}{j}"
        strikes = strikes_r[trialname] if side == "right" else strikes_l[trialname]
        ax = figures[str(i)].add_subplot(3, int(numtrials/3), j)
        angle_plot(trialname, strikes, label, ax)
        subplotsdict[trialname] = ax

    figures[str(i)].suptitle(f"{trialtitle} {side.capitalize()} {jointlabel}")

    # Set consistent y-limits
    allYLim = np.array([ax.get_ylim() for ax in subplotsdict.values()])
    minYLim, maxYLim = allYLim.min(), allYLim.max()
    for ax in subplotsdict.values():
        ax.set_ylim(minYLim, maxYLim)

    figures[str(i)].tight_layout()

# Show all figures
plt.show()