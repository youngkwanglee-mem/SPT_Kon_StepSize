% SPT_Kon_StepSize.m
% Version 1.0
% This script loads TrackMate CSV output files to analyze association 
% kinetics and create step size distributions.
% The folder must contain one raw TIF image analyzed by TrackMate
% and a 'spots.csv' file exported from TrackMate. Subfolders containing 
% additional 'spots.csv' files are not permitted.
% For TrackMate analysis, use the default settings for units (1 pixel) 
% and time (1 frame).
% 
% MATLAB version: 2022b Update 5
% TrackMate version: 7.10.2

%% FILE PATH
clc
clear all

% Define the raw TIF image file path. The directory must contain 
% exactly one image and one 'spots.csv' file.
imagefile = '/Users/YOUNGKWANG/Desktop/code cleaning/06_S12_2pct_PS_Atto488_DPPE_9_positions_369_frames.tif'
[rootdir, imagefilename, imageext] = fileparts(imagefile);

filelist = dir(fullfile(rootdir, 'spots.csv'));
InfoImage = imfinfo(imagefile);
filelist = filelist(~[filelist.isdir]); % Remove directories from the list

for k=1:length(filelist);
    spots_csv_file{k} = fullfile(filelist(k).folder, filelist(k).name);
end

NumFile = size(spots_csv_file,2);
l = 1; % Single-file analysis mode

%% USER INPUTS
pixelsize = 0.108; % um
dt = 0.1085; % sec
max_displacement_trackmate = 10; % pixel; The maximum displacement used during TrackMate analysis.
num_positions = 1; % Number of concatenated positions. Enter 1 for a single position image stack.

% Trajectory filtering parameters:
% Mean speed filter
mean_velocity_low_threshold = 0.5; % pixel per frame; identical to the "track mean speed" threshold in TrackMate.
diff_min_num_frames = 5; % frames; Minimum track length for step size distribution analysis. 2 includes all links.

% Option to remove transiently immobile tracks
window = 3; % Number of consecutive immobile steps required for removal.
cutoff_displacement = 0; % pixel; Enter 0 to disable. Otherwise, enter the minimum displacement for the window steps.

% Option to remove overlapping particles (colliding tracks)
remove_colliding_track = 0; % Enter 0 to disable. Enter 1 to remove trajectories that approach closer than the 'collision_cutoff_distance'.
collision_cutoff_distance = 6; % pixel; Distance threshold for identifying colliding spots.

% ROI settings
roi_width = 300; % Square ROI width. Should be at least 10 pixels smaller than full width to exclude edge particles.
roi_shift_x = 0; % Shift from center (0 = centered)
roi_shift_y = 0; % Shift from center (0 = centered)
full_width = InfoImage(1).Width; % pixel

% Other analysis parameters
stepsize_binning_um = 0.01; % um; Bin size for the step size distribution (x-axis).

%% PARAMETER CALCULATION AND LOGGING
num_frame = length(InfoImage); % Total number of frames in the image
num_frame_per_pos = num_frame/num_positions; % Frames per position

% Define ROI boundaries
position_xleft_cutoff = full_width/2+roi_shift_x-roi_width/2-1;
position_xright_cutoff = full_width/2+roi_shift_x+roi_width/2-1;
position_ytop_cutoff = full_width/2+roi_shift_y-roi_width/2-1;
position_ybottom_cutoff = full_width/2+roi_shift_y+roi_width/2-1;

% Convert max displacement to microns
links_maxdistance = max_displacement_trackmate*pixelsize; 

% Organize analysis parameters into a table for export
Analysis_parameters = [pixelsize max_displacement_trackmate dt num_frame num_positions num_frame_per_pos window cutoff_displacement diff_min_num_frames...
    position_xleft_cutoff position_xright_cutoff position_ytop_cutoff position_ybottom_cutoff...
    stepsize_binning_um roi_width roi_shift_x roi_shift_y remove_colliding_track collision_cutoff_distance mean_velocity_low_threshold];

T_Analysis_parameters = array2table(Analysis_parameters);
T_Analysis_parameters.Properties.VariableNames(1:size(Analysis_parameters,2)) = ...
    {'pixel size (um)', 'max displacement in trackmate', 'dt (sec)', 'total number of frames', ...
    'number of positions', 'number of frames per position', 'stationary window steps', ...
    'displacement cutoff for stationary window', 'minimum frames for step size distribution', ...
    'position_xleft_cutoff', 'position_xright_cutoff', 'position_ytop_cutoff', 'position_ybottom_cutoff', ...
    'step size bin(um)', 'roi_width', 'roi_shift_x', 'roishift_y', ...
    'remove colliding track (1:yes, others:no)', 'colliding threshold (pixel)', 'Mean velocity threshold (pixel)'};

writetable(rows2vars(T_Analysis_parameters), fullfile(filelist(l).folder,'output_00_Analysis_parameters.xlsx'));

start_frames_per_pos = [0:num_frame_per_pos:num_frame-1]; % Frame indexing starts at 0

%% LOAD DATA
file_path_spots_csv = fullfile(filelist(l).folder, filelist(l).name);
spot_csv_org = readtable(file_path_spots_csv,'ReadRowNames',true,...
    'VariableDescriptionsLine',2,'VariableUnitsLine',4);
spot_csv_org = spot_csv_org(2:end,1:end); % Remove the first row containing N/A metadata

%% REMOVE EDGE TRACKS
% Extract single-frame detections
spot_csv_single_idx = isnan(spot_csv_org.TRACK_ID);
spot_csv_single = spot_csv_org(spot_csv_single_idx,:);

% Extract multi-frame track data
spot_csv_multi_idx = ~isnan(spot_csv_org.TRACK_ID);
spot_csv_multi = spot_csv_org(spot_csv_multi_idx,:);

% Calculate the mean X and Y coordinates for each multi-frame track
track_id_multi = unique(spot_csv_multi.TRACK_ID);
track_meanXY = zeros(size(track_id_multi,1),3);
for n = 1:size(track_id_multi)
    idx = find(spot_csv_multi.TRACK_ID == track_id_multi(n));
    meanX = mean(spot_csv_multi.POSITION_X(idx));
    meanY = mean(spot_csv_multi.POSITION_Y(idx));
    track_meanXY(n,:) = [track_id_multi(n) meanX meanY];
end

tbl_track_meanXY = table(track_meanXY(:,1), track_meanXY(:,2), track_meanXY(:,3));
tbl_track_meanXY.Properties.VariableNames = ["TRACK_ID", "MEAN_X", "MEAN_Y"];

% Filter multi-frame particles outside the ROI
edge_track_idx = tbl_track_meanXY.MEAN_X < position_xleft_cutoff | tbl_track_meanXY.MEAN_X > position_xright_cutoff | tbl_track_meanXY.MEAN_Y < position_ytop_cutoff | tbl_track_meanXY.MEAN_Y > position_ybottom_cutoff;
edge_track_ID = unique(tbl_track_meanXY.TRACK_ID(edge_track_idx));
center_track_ID = unique(tbl_track_meanXY.TRACK_ID(~edge_track_idx));
center_trackID_idx = ismember(spot_csv_multi.TRACK_ID, center_track_ID);

spot_csv_multi_edge_removed = spot_csv_multi(center_trackID_idx,:); 
spot_csv_multi_edge_removed = sortrows(spot_csv_multi_edge_removed,{'TRACK_ID','FRAME'});

% Filter single-frame particles outside the ROI
edge_spot_idx = spot_csv_single.POSITION_X < position_xleft_cutoff | spot_csv_single.POSITION_X > position_xright_cutoff | spot_csv_single.POSITION_Y < position_ytop_cutoff | spot_csv_single.POSITION_Y > position_ybottom_cutoff;
spot_csv_single_edge_remove = spot_csv_single(~edge_spot_idx,:);

% Combine filtered single and multi-frame data
spot_csv_edge_removed = [spot_csv_multi_edge_removed; spot_csv_single_edge_remove];
center_track_ID = unique(spot_csv_edge_removed.TRACK_ID);

%% REMOVE COLLIDING (NEIGHBORING) TRAJECTORIES
if remove_colliding_track == 1;
    % Identify tracks that approach each other too closely
    first_frame = min(spot_csv_multi_edge_removed.FRAME);
    last_frame = max(spot_csv_multi_edge_removed.FRAME);
    colliding_track_ID = double.empty();
    
    for q = first_frame:last_frame
        FRAME_q_idx = find(spot_csv_multi_edge_removed.FRAME == q);
        spots_q = spot_csv_multi_edge_removed(FRAME_q_idx,:);
        
        if size(spots_q,1) > 1
            spot_pair_idx = nchoosek(1:size(spots_q,1),2); % Create pairing indices for coordinates
            x_pairs = spots_q.POSITION_X(spot_pair_idx);
            
            % Handle dimension mismatch when only two particles exist
            if size(x_pairs) == [2 1]
                x_pairs = x_pairs';
            end
            
            pair_distance(:,1) = diff(x_pairs,1,2); % Pixel distance between pairs (X)
            
            y_pairs = spots_q.POSITION_Y(spot_pair_idx);
            if size(y_pairs) == [2 1]
                y_pairs = y_pairs';
            end
            
            pair_distance(:,2) = diff(y_pairs,1,2); % Pixel distance between pairs (Y)
            pair_distance = abs(pair_distance);
            
            % Find indices of neighbors based on the collision threshold
            x_neighbours = find(pair_distance(:,1) < collision_cutoff_distance);
            y_neighbours = find(pair_distance(:,2) < collision_cutoff_distance);
            true_colliding_spot_pair_idx = intersect(x_neighbours, y_neighbours);
            
            colliding_spot_idx = unique(spot_pair_idx(true_colliding_spot_pair_idx,:));
            colliding_track_ID_to_add = unique(spots_q.TRACK_ID(colliding_spot_idx));
            colliding_track_ID = [colliding_track_ID; colliding_track_ID_to_add];
            
            pair_distance = double.empty(); % Reset variable to prevent dimension errors
        end
    end
    
    colliding_track_ID = unique(colliding_track_ID);
    num_colloiding_track = size(colliding_track_ID,1);
    fprintf('%d colliding tracks were removed.\n', num_colloiding_track)
    
    % Filter out colliding tracks from the table
    isolated_trajectory_idx = ~ismember(spot_csv_multi_edge_removed.TRACK_ID, colliding_track_ID);
    spot_csv_multi_edge_colliding_removed = spot_csv_multi_edge_removed(isolated_trajectory_idx,:);
    center_isolated_track_ID = unique(spot_csv_multi_edge_colliding_removed.TRACK_ID);
else
    spot_csv_multi_edge_colliding_removed = spot_csv_multi_edge_removed;
    center_isolated_track_ID = center_track_ID;
    colliding_track_ID = double.empty();
    fprintf('Collision removal was disabled.\n')
end

%% REMOVE IMMOBILE AND SLOW TRAJECTORIES
% Calculate displacements
displacement_x = diff(spot_csv_multi_edge_colliding_removed.POSITION_X);
displacement_y = diff(spot_csv_multi_edge_colliding_removed.POSITION_Y);
displacement_r = sqrt(displacement_x.^2 + displacement_y.^2); 

% Flag track starts to avoid calculating displacement between different tracks
track_start = diff(spot_csv_multi_edge_colliding_removed.TRACK_ID);
track_start_idx = find(track_start ~= 0);
displacement_r(track_start_idx) = nan; 
displacement_r = [nan; displacement_r]; % Pad to match original matrix dimensions
spot_csv_multi_edge_colliding_removed.DISPLACEMENT = displacement_r;

% Detect immobile displacements using a sliding window
for j=1:size(spot_csv_multi_edge_colliding_removed,1)-window+1;
    filter_logic_test = spot_csv_multi_edge_colliding_removed.DISPLACEMENT(j:j+window-1) < cutoff_displacement;
    filter_logic_test(window+1) = 1; % Dummy value to prevent errors on empty sets
    counts = sum(filter_logic_test);
    filter_track_idx(j,1) = (counts == window+1);
end

% Identify trajectories containing immobile segments
immobile_trajectory_ID = spot_csv_multi_edge_colliding_removed.TRACK_ID(1:size(filter_track_idx,1),1);
immobile_trajectory_ID = unique(immobile_trajectory_ID(filter_track_idx));

total_tracks_number = size(unique(spot_csv_multi_edge_colliding_removed.TRACK_ID),1);
immobile_trajectory_number = size(immobile_trajectory_ID,1);
fprintf('%d immobile tracks were detected.\n', immobile_trajectory_number);

% Filter immobile tracks from the table
mobile_trajectory_idx = ~ismember(spot_csv_multi_edge_colliding_removed.TRACK_ID, immobile_trajectory_ID);
spot_csv_multi_edge_colliding_immobile_removed = spot_csv_multi_edge_colliding_removed(mobile_trajectory_idx,:);
center_isolated_mobile_track_ID = unique(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID);

% Identify tracks with mean velocity below the threshold (slow tracks)
% This subset is specifically for diffusion analysis
mobile_track_ID = unique(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID);
mean_velocity = zeros(size(mobile_track_ID));
for s = 1: size(mobile_track_ID);
    idx = find(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID == mobile_track_ID(s));
    displacement = spot_csv_multi_edge_colliding_immobile_removed.DISPLACEMENT(idx);
    mean_velocity(s) = mean(displacement(~isnan(displacement))); 
end

not_slow_track_ID = mobile_track_ID(mean_velocity > mean_velocity_low_threshold);
slow_track_ID = mobile_track_ID(mean_velocity < mean_velocity_low_threshold);
fprintf('%d slow tracks were detected.\n', size(slow_track_ID,1));

% Final table for normal diffusion analysis (mobile and not slow)
mobile_notSlow_trajectory_idx = ismember(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID, not_slow_track_ID);
spot_csv_multi_edge_colliding_immobile_slow_removed = spot_csv_multi_edge_colliding_immobile_removed(mobile_notSlow_trajectory_idx,:);
center_isolated_mobile_notSlow_track_ID = unique(spot_csv_multi_edge_colliding_immobile_slow_removed.TRACK_ID);

%% RECRUITMENT ANALYSIS
% Calculate the recruitment rate by identifying the first appearance of each track
mult_dwell_first_frm = zeros(size(center_isolated_mobile_track_ID,1),1);
mult_dwell_first_frm_trackid = zeros(size(center_isolated_mobile_track_ID,1),1);

for n = 1:size(center_isolated_mobile_track_ID)
    mult_dwell_first_frm(n) = min(spot_csv_multi_edge_colliding_immobile_removed.FRAME(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID == center_isolated_mobile_track_ID(n)));
    mult_dwell_first_frm_trackid(n) = center_isolated_mobile_track_ID(n);
end

% Remove pre-bound particles (those appearing in the very first frame of a stack)
mult_dwell_first_frm = mult_dwell_first_frm(~ismember(mult_dwell_first_frm, start_frames_per_pos));
single_dwell_first_frm = spot_csv_single_edge_remove.FRAME(~ismember(spot_csv_single_edge_remove.FRAME, start_frames_per_pos));

% Merge single-frame and multi-frame recruitment events
first_frm = [single_dwell_first_frm; mult_dwell_first_frm];

% Calculate the cumulative histogram of recruitment events
recruitment_time_frame = 1:num_frame;
cumulative_recruitment_all = cumsum(histc(first_frm, recruitment_time_frame));
recruitment_time_sec = recruitment_time_frame.*dt - dt/2; % Calculate bin centers

% Fit a linear model (intercept forced to zero) to determine the recruitment rate
rate = fitlm(recruitment_time_sec, cumulative_recruitment_all, 'Intercept', false);
rate_s_um2 = table2array(rate.Coefficients(1,1)) / (roi_width*pixelsize)^2;

% Repeat recruitment calculation for multi-dwell particles only
cumulative_recruitment_multidwell_only = cumsum(histc(mult_dwell_first_frm, recruitment_time_frame));
rate_multi = fitlm(recruitment_time_sec, cumulative_recruitment_multidwell_only, 'Intercept', false);
rate_s_um2_multi = table2array(rate_multi.Coefficients(1,1)) / (roi_width*pixelsize)^2;

%% SAVE RECRUITMENT DATA
recruitment_plot_save = array2table([recruitment_time_sec' cumulative_recruitment_all rate.Fitted cumulative_recruitment_multidwell_only rate_multi.Fitted ],...
    'VariableNames', {'time_sec', 'cumulative_recruitment_per_sec', 'fit', 'cumulative_recruitment_no_single_dwell', 'fit_no_single_dwell'});

recruitment_fit_coefficients = rate.Coefficients;
recruitment_fit_coefficients.Rate = rate_s_um2;
recruitment_fit_coefficients.Properties.VariableNames("Estimate") = "Frequency_s_inv";
recruitment_fit_coefficients.Properties.VariableNames("Rate") = "Rate_s_inv_um_sq_inv";

writetable(recruitment_plot_save, fullfile(filelist(l).folder,'output_01-1_Recruitment_data_fit.xlsx'));
writetable(recruitment_fit_coefficients, fullfile(filelist(l).folder,'output_01-2_Recruitment_fit_coefficients.xlsx'));

%% STEP SIZE DISTRIBUTION FOR DIFFUSION ANALYSIS
% Exclude pre-bound particles from dwell time calculations
dwell_TRACK_ID = mult_dwell_first_frm_trackid(~ismember(mult_dwell_first_frm, start_frames_per_pos));
dwell_frame_length_multi = zeros(size(dwell_TRACK_ID,1),1);

for m = 1:size(dwell_TRACK_ID,1);
    dwell_TRACK_ID_logic = ismember(spot_csv_multi_edge_colliding_immobile_removed.TRACK_ID, dwell_TRACK_ID(m));
    dwell_frames = spot_csv_multi_edge_colliding_immobile_removed.FRAME(dwell_TRACK_ID_logic);
    dwell_frame_length_multi(m) = max(dwell_frames) - min(dwell_frames) + 1;
end

% Filter out tracks shorter than the specified minimum frame count
track_ID_step_size_logic = dwell_frame_length_multi >= diff_min_num_frames;
track_ID_step_size = dwell_TRACK_ID(track_ID_step_size_logic);

% Collect displacements for step size analysis and convert to microns
step_size_um = spot_csv_multi_edge_colliding_immobile_slow_removed.DISPLACEMENT(ismember(spot_csv_multi_edge_colliding_immobile_slow_removed.TRACK_ID, track_ID_step_size)).*pixelsize;
step_size_um = step_size_um(~isnan(step_size_um));

% Generate step size histogram and probability density
step_size_um_bin = [0.0001:stepsize_binning_um:max(step_size_um)+0.0001]'; 
step_size_um_freq = histc(step_size_um, step_size_um_bin);
step_size_um_probd = step_size_um_freq / sum(step_size_um_freq) / stepsize_binning_um;

%% SAVE STEP SIZE DISTRIBUTION
step_size_save = array2table([step_size_um_bin step_size_um_probd],...
    'VariableNames', {'step_size_um', 'probability_density'});
writetable(step_size_save, fullfile(filelist(l).folder,'output_02_Step_size_data.xlsx'));

%% SAVE FILTERED SPOTS AND BACKUP CODE
writetable(spot_csv_multi_edge_colliding_immobile_slow_removed, fullfile(filelist(l).folder,'output_03_Filtered_spots.xlsx'));
save(fullfile(rootdir,'output_00_All_variables.mat'));

% Create a backup of the current script in the output directory
code_file = [mfilename('fullpath')];
[codefile_path, codefile_name, codefile_ext] = fileparts(code_file);
newbackup = sprintf('%s.m', convertCharsToStrings(rootdir) + "/output_00_CODE_BACKUP_" + convertCharsToStrings(codefile_name));
currentfile = strcat(code_file, '.m');
copyfile(currentfile, newbackup);

disp('Analysis completed.')