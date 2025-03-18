%clear;

%% select files
% assumptions:
%   the two files are not synchronized
%   joystick_data_filename contains raw joystick data and event markers
%   spikes_data_filename contains spikes (perhaps sorted) and event markers

joystick_data_filename = 'test_data\Q_20201109_ICMS_dim5.plx';
spikes_data_filename = 'test_data\Q_20201109_ICMS_dim5.nev';
spikes_channel = 95;
spikes_unit = 3;

%% select patterns
% the basic pattern is for any successful COT trial

% start trial marker block
% 	// 1: EventCode::TrialStart = 251
% 	// 2: serial 0
% 	// 3: serial 1
% 	// 4: serial 2
% 	// 5: serial 3
% 	// 6: input mode
% 	// 7: trial kind
% 	// 8: annuli count
% 	// 9: targets count
% 	//10: annulus index
% 	//11: target index
% 	//12: prior success (EventCode::PriorSuccess or EventCode::PriorFailure)
% 	//13: trial execution mode (EventCode::ExecuctionMode or EventCode::ObservationMode)

start_trial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.');
basic_pattern = Pattern(start_trial, 210, 245, 60, 248, 119, 120, 121, 250, 253);

% 4 patterns with 4 targets in cardinal directions
patterns = [basic_pattern.replace_at(1, start_trial.replace_at(11, 0)),...
    basic_pattern.replace_at(1, start_trial.replace_at(11, 2)),...
    basic_pattern.replace_at(1, start_trial.replace_at(11, 4)),...
    basic_pattern.replace_at(1, start_trial.replace_at(11, 6))];

%patterns = Pattern(start_trial, 210, 245, 60, 248, 119, 120, 121, 250, 253);

%% joystick processing

js_sp = FMLReadSpeedTimeseries(joystick_data_filename);

% read marker events
[js_markers_ts, js_markers_sv] = FMLReadMarkers(joystick_data_filename);

% extract trials and align on first marker
[js_T, js_MT, js_patsIdx, js_patsT ] = TimeseriesFindExtractAlign(js_sp, js_markers_ts, js_markers_sv, patterns);

% calculate movement onset
js_onsets = FMLFindOnset(js_T, js_MT);

% to modify FMLFindOnset options:
% onset_options = FMLFindOnset; % default options
% onset_options.baselinefromindex = 2;  % index of marker that anchors the baseline calculation
% onset_options.baselineduration = [-0.5, 0]; % from 500ms before baseline marker time, to 0ms before the baseline marker
% onset_options.onsetfromindex = 3;    % index of the marker from which time we seach for onset
% onset_options.onsetstddevfactor = 3; % onset is when over this factor times stddev of baseline from baseline...
% onset_options.onsetduration = 0.03; % ... from this amount of time
% onset_options.reachtoindex = 5; % to calculate movement duration: from onset to this marker
% js_onsets = FMLFindOnset(js_T, js_MT, onset_options);

% optionally plot some histograms
plot_onsets_histograms = true;
% plot_onsets_histograms = true;
if plot_onsets_histograms
    for i=1:length(js_onsets)
        figure;
        histogram([js_onsets{i}(:).reaction_duration], 100);
        title(['pattern(' num2str(i) ') reaction duration']);
        figure;
        histogram([js_onsets{i}(:).movement_duration], 100);
        title(['pattern(' num2str(i) ') movement duration']);
    end
end

%% spikes processing

% read spikes and markers
% the reason for the various checks (e.g., ~exist, etc) is to prevent
% re-reading the same file. FMLReadSpikeTimes and FMLReadMarkers can be
% slow to read on .nev files, for example.
if ~exist('sp_spikes_ts', 'Var') || ~exist('old_spikes_data_filename', 'Var') || ...
        ~strcmp(spikes_data_filename, old_spikes_data_filename) || ~exist('old_spikes_channel', 'Var') || ...
        ~exist('old_spikes_unit', 'Var') || old_spikes_channel ~= spikes_channel || old_spikes_unit ~= spikes_unit
    old_spikes_data_filename = spikes_data_filename;
    old_spikes_channel = spikes_channel;
    old_spikes_unit = spikes_unit;
    [sp_spikes_ts] = FMLReadSpikeTimes(spikes_data_filename, spikes_channel, spikes_unit);
    [sp_markers_ts, sp_markers_sv] = FMLReadMarkers(spikes_data_filename);
end

% extract trials and align on first marker
if ~isempty(sp_spikes_ts)
    [sp_T, sp_MT, sp_patsIdx, sp_patsT ] = SpikesFindExtractAlign(sp_spikes_ts, sp_markers_ts, sp_markers_sv, patterns);
else
    disp('empty sp_spikes_ts');
end

%% plot raster histograms aligned on movement onset

% to re-align on movement onset we shift by the onset time 
% also insert a marker at 0 (movement onset marker)
js_T_r = cell(1,length(js_T));
js_MT_r = cell(1,length(js_T));
sp_T_r = cell(1,length(js_T));
sp_MT_r = cell(1,length(js_T));

for i=1:length(sp_T)
    for it=1:length(sp_T{i})
        js_T_r{i}{it}.StartTime = js_T{i}{it}.StartTime - js_onsets{i}(it).onset_time;
        js_MT_r{i}{it} = sort([js_MT{i}{it} - js_onsets{i}(it).onset_time, 0]);
        sp_T_r{i}{it} = sp_T{i}{it} - js_onsets{i}(it).onset_time;
        sp_MT_r{i}{it} = sort([sp_MT{i}{it} - js_onsets{i}(it).onset_time, 0]); 
    end
end

plot_rasters = true;
if plot_rasters
    figure;
    SpikesMultiRasterHistogram(sp_T_r, sp_MT_r, 100);
end

