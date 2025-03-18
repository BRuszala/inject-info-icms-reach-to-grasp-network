%% Add paths to necessary toolboxes

clear; %close all; 
addpath(genpath('F:\Projects\COT_ICMS\distantMod\'))
addpath(genpath('F:\MATLAB\FML_matlab'))
addpath(genpath('F:\MATLAB\fileHandling_BR'))

%% Call Files to analyze

monk ='Qulio';
cortical_region = ['PMd'];
useTitle = false;
useResT = false;
useCatchT = false;
ResT = 3;
CatchT = 3;

path_to_data = ['F:\RawData\COT_ICMS\' monk '\' cortical_region '\'];
savePath = ['F:\Publications\InjectInformationPMCPPC\Results\'];
all_nevFiles = [dir([path_to_data 'Ripple\' monk(1) '_*ICMS_*.nev'])];% dir([path_to_data 'Ripple\' monk(1) '_*ICMS_baseline*.nev'])];
all_ns2Files = [dir([path_to_data 'Ripple\' monk(1) '_*ICMS_*.ns2'])];% dir([path_to_data 'Ripple\' monk(1) '_*ICMS_baseline*.ns2'])];
all_h5Files = [dir([path_to_data 'UTCS\' monk(1) '_*ICMS_*.h5'])];% dir([path_to_data 'UTCS\' monk(1) '_*ICMS_baseline*.h5'])];
nev_fileList = {all_nevFiles(:).name}';
ns2_fileList = {all_ns2Files(:).name}';
h5_fileList = {all_h5Files(:).name}';

% Check for existing .mat file with previous information and list of used files
% If there is no .mat file a new one will be created

if exist([path_to_data monk '_RT_MT_Success_' cortical_region '.mat'], 'file') ~= 0 
    load([path_to_data monk '_RT_MT_Success_' cortical_region '.mat'], 'rt_mt_perform_Info')           
    for fi = 1:length(rt_mt_perform_Info)
        try
            [thisPath, thisFile, thisExt] = fileparts(rt_mt_perform_Info(fi).File);           
            usedFiles{fi,1} = [thisFile, thisExt];
        catch
            usedFiles{fi,1} = [];
        end
    end
    newNDX = cellfun(@(x) ~sum(strcmpi(x, usedFiles)), ns2_fileList);     
else    
    rt_mt_perform_Info = struct('File', [], 'Date', [], 'targetRT', [], 'targetMT', [], 'targetSuccess', [],...
                            'targetTrials', []);
    newNDX = true(length(ns2_fileList),1);    
end

                        
%% Build patterns
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

start_trial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 31, '.');  %start_trial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.');
basic_pattern = Pattern(start_trial, 210, 245, 60, 248, 119, 120, 121, 250, 253);

% start_trial = Marker(251, '.', '.', '.', '.', '.', 10, '.', '.', '.', '.', 31, '.');  %start_trial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.');
% basic_pattern = Pattern(start_trial, 210, 245, 60, 253);


% 4 patterns with 4 targets in cardinal directions
patterns8 = [basic_pattern.replace_at(1, start_trial.replace_at(11, 0)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 2)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 4)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 6))];
patterns4 = [basic_pattern.replace_at(1, start_trial.replace_at(11, 0)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 1)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 2)),...
            basic_pattern.replace_at(1, start_trial.replace_at(11, 3))];
        
if useCatchT % Treat catch trials like normal trials
    catchPattern = basic_pattern.replace_at(1, start_trial.replace_at(11, catchT));
    patterns8 = [patterns8, catchPattern];
    patterns4 = [patterns4, catchPattern];
end     
if useResT % A "success" pattern here simply means the monkey left the center (for electrode resolution studies)
    ResPattern = Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', ResT, 31, '.'), 210, 245, 60, 253);
    patterns8 = [patterns8, ResPattern];
    patterns4 = [patterns4, ResPattern];
end       
successPatterns8 = patterns8;
successPatterns4 = patterns4;

trialPatterns8 = [Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 0, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 2, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 4, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 6, 31, '.'), Marker(253))];

trialPatterns4 = [Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 0, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 1, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 2, 31, '.'), Marker(253)),...
                 Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', 3, 31, '.'), Marker(253))];
             
if useCatchT
    catchTrial = Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', catchT, 31, '.'), Marker(253));
    trialPatterns8 = [trialPatterns8, catchTrial];
    trialPatterns4 = [trialPatterns4, catchTrial];
end              
if useResT
    resTrial = Pattern(Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', ResT, 31, '.'), Marker(253));
    trialPatterns8 = [trialPatterns8, resTrial];
    trialPatterns4 = [trialPatterns4, resTrial];
end 
%patterns = Pattern(start_trial, 210, 245, 60, 248, 119, 120, 121, 250, 253);


%% Loop through list of new files

% newNDX(:) = true;
%[ns2_fileList, nev_fileList, [h5_fileList; {''}]]
[ns2_fileList, nev_fileList, h5_fileList]
nev_fileList = nev_fileList(newNDX);
ns2_fileList = ns2_fileList(newNDX);
h5_fileList = h5_fileList(newNDX);
[ns2_fileList, nev_fileList, h5_fileList]
correctEventCodes = false(size(ns2_fileList));
numUsedFiles = sum(~cellfun(@(x) isempty(x), {rt_mt_perform_Info(:).File}'));
%%

startTime = tic;
for fi = 1:length(ns2_fileList)
    
    disp(['File ' num2str(fi) ' of ' num2str(sum(newNDX)) ' new files'])
    rt_mt_NDX = fi+numUsedFiles;
    joystick_data_filename = [path_to_data 'Ripple\' ns2_fileList{fi}];
    marker_data_filename = [path_to_data 'Ripple\' nev_fileList{fi}];
    h5_filename = [path_to_data 'UTCS\' h5_fileList{fi}];
    
    parsedFilename = regexp(ns2_fileList{fi}, '_', 'split');
    date = parsedFilename{2};

    rt_mt_perform_Info(rt_mt_NDX).File = joystick_data_filename;
    rt_mt_perform_Info(rt_mt_NDX).Date = date;

    if (strcmpi(monk, 'Felix') && str2num(date) >= 20220104) || strcmpi(monk, 'Indira')
        trialPatterns = trialPatterns4;
        successPatterns = successPatterns4;
        patterns = patterns4;
    else
        trialPatterns = trialPatterns8;
        successPatterns = successPatterns8;
        patterns = patterns8;
    end
    
    %%%%% Joystick processing %%%%%

    try
        js_sp = FMLReadSpeedTimeseries(joystick_data_filename, 'analog 10','analog 11');
    catch
        warning('Using Analog 3&4')
        js_sp = FMLReadSpeedTimeseries(joystick_data_filename, 'analog 3','analog 4');
    end

    % read marker events
    [js_markers_ts, js_markers_sv] = FMLReadMarkers(marker_data_filename);

    if correctEventCodes(fi)

        clear trial_data
        %Reading H5 file
        h5_filename = [path_to_data 'UTCS\' h5_fileList{fi}];
        h5_file_info = h5info(h5_filename);
        for trial = 1:length(h5_file_info.Groups)
            trial_info = h5info(h5_filename, h5_file_info.Groups(trial).Name);

            % Extract information saved in trial Attributes
            for attribute_ind = 1:size(trial_info.Attributes,1)
                trial_data(trial).(genvarname(trial_info.Attributes(attribute_ind).Name)) = ...
                    trial_info.Attributes(attribute_ind).Value;
            end

            % load trial information under root trial group
            for dataset_ind = 1:length(trial_info.Datasets)
                trial_data(trial).(genvarname(trial_info.Datasets(dataset_ind).Name)) = ...
                    h5read(h5_filename, [h5_file_info.Groups(trial).Name '/' trial_info.Datasets(dataset_ind).Name]);
                % generate variable name using trial_info.Datasets(dataset_ind).Name
                % trial_data stores all datasets from BCI controlled test
            end

        end

        marker_data = js_markers_sv;
        marker_ts = js_markers_ts;
        trueCodes = cat(1,trial_data.event_codes);
        [filterCodes, filterTS] = fixEventCodes_BR(trueCodes, marker_data, marker_ts);

        js_markers_sv = filterCodes;
        js_markers_ts = filterTS;

    end

    try
        trials_count = h5readatt(h5_filename, '/', 'trials_count');
        ec = cell(trials_count,1);
        blockNames = string(nan(size(ec)));
        h5_file_info = h5info(h5_filename);
        for i=1:trials_count
            ec{i} = h5read(h5_filename, ['/trial_' pad(num2str(i), 4, 'left', '0'), '/event_codes']);

            trial_info = h5info(h5_filename, h5_file_info.Groups(i).Name);
            blockNames(i) = string(trial_info.Attributes(8).Value);
        end
        [clrs, changes]= COTReportTargetColorChange(h5_filename);
        clrs = translateClrs(clrs);
%         clrs(strcmpi(blockNames, "5e_noStim")) = "dimFull";
        theseClrs = unique(clrs);
        existClrsAtt = true;
    catch
        existClrsAtt = false;
        warning('No Target Color Attribute')
        theseClrs = "No Color Attribute in .h5";
    end

    % extract trials and align on first marker
    [js_T, js_MT, js_patsIdx, js_patsT] = TimeseriesFindExtractAlign(js_sp, js_markers_ts, js_markers_sv, successPatterns);

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

    %%%%% Success processing %%%%%

    successPerTarget = FindPattern(js_markers_sv, successPatterns);
    trialsPerTarget = FindPattern(js_markers_sv, trialPatterns);

    successID = cellfun(@(x) getTrialID(x,js_markers_sv), successPerTarget, 'Uni', 0);
    trialID = cellfun(@(x) getTrialID(x,js_markers_sv), trialsPerTarget, 'Uni', 0);

    freqClrs = zeros(size(theseClrs));
    trialsPerTargetColor = zeros(length(theseClrs), length(trialsPerTarget));
    successPerTargetColor = zeros(size(trialsPerTargetColor));
    rtPerTargetColor = cell(size(trialsPerTargetColor));
    mtPerTargetColor = cell(size(trialsPerTargetColor));
    
    if existClrsAtt                 
        for targ_i = 1:length(successPerTarget)           
            successClrs = clrs(cellfun(@(x) any(sum(x(find(x==251,1,'first')+1:find(x==251,1,'first')+4)' == successID{targ_i}, 2) == 4), ec));
            trialsClrs = clrs(cellfun(@(x) any(sum(x(find(x==251,1,'first')+1:find(x==251,1,'first')+4)' == trialID{targ_i}, 2) == 4), ec));
            if ~isempty(trialsClrs)
                for ci = 1:length(theseClrs)

                    trialColorNDX = strcmpi(trialsClrs, theseClrs(ci));
                    successColorNDX = strcmpi(successClrs, theseClrs(ci));
                    freqClrs(ci) = sum(strcmpi(trialsClrs, theseClrs(ci)))/length(trialsClrs);
                    trialsPerTargetColor(ci, targ_i) = size(trialsPerTarget{targ_i}(trialColorNDX, :),1);
                    successPerTargetColor(ci, targ_i) = size(successPerTarget{targ_i}(successColorNDX, :),1);
                    rtPerTargetColor{ci, targ_i} = [js_onsets{targ_i}(successColorNDX).reaction_duration]';
                    mtPerTargetColor{ci, targ_i} = [js_onsets{targ_i}(successColorNDX).movement_duration]';

                end
            else
                trialsPerTargetColor(1, targ_i) = NaN;
                successPerTargetColor(1, targ_i) = NaN;
                rtPerTargetColor{1, targ_i} = NaN;
                mtPerTargetColor{1, targ_i} = NaN;
            end
        end
    else
        freqClrs = "NaN";
        for targ_i = 1:length(successPerTarget)
            trialsPerTargetColor(1, targ_i) = size(trialsPerTarget{targ_i}(:, :),1);
            successPerTargetColor(1, targ_i) = size(successPerTarget{targ_i}(:, :),1);
            rtPerTargetColor{1, targ_i} = [js_onsets{targ_i}(:).reaction_duration]';
            mtPerTargetColor{1, targ_i} = [js_onsets{targ_i}(:).movement_duration]';
        end
    end
    
    rt_mt_perform_Info(rt_mt_NDX).targetTrials = trialsPerTargetColor;
    rt_mt_perform_Info(rt_mt_NDX).targetSuccess = successPerTargetColor;
    rt_mt_perform_Info(rt_mt_NDX).targetRT = rtPerTargetColor;
    rt_mt_perform_Info(rt_mt_NDX).targetMT = mtPerTargetColor;
    rt_mt_perform_Info(rt_mt_NDX).allColors = theseClrs';
    rt_mt_perform_Info(rt_mt_NDX).colorFreqs = freqClrs';

end
endTime = toc(startTime);
%%
for fi = 1:length(rt_mt_perform_Info)   
    thisFile = rt_mt_perform_Info(fi).File;
    if ~isempty(thisFile)        
        slashNDX = regexp(thisFile, '\');
        scoreNDX = regexp(thisFile, '_');
        dateNDX = scoreNDX(scoreNDX > slashNDX(end));
        allDates(fi,:) = string(thisFile(dateNDX(1)+1:dateNDX(2)-1));
    else
        continue
    end
end

%Add new files into existing list crhonilogically by date

[~, sortNDX] = sort(allDates);
rt_mt_perform_Info = rt_mt_perform_Info(sortNDX);

% load(['F:\Projects\COT_ICMS\catchTrialRts.mat'], 'catchInfo')
% catchSession = 46;
% catchInfo(end+1).monk = monk;
% catchInfo(end).corticalRegion = cortical_region;
% catchInfo(end).rts = vertcat(rt_mt_perform_Info(catchSession).targetRT{:});
% save(['F:\Projects\COT_ICMS\catchTrialRts.mat'], 'catchInfo')

save([path_to_data monk '_RT_MT_Success_' cortical_region '.mat'], 'rt_mt_perform_Info') %Save information


%% Plotting Data
close all

% monk = 'Qulio';
% region = 'AIP';
path_to_data = ['F:\RawData\COT_ICMS\' monk '\' cortical_region '\'];

load([path_to_data monk '_RT_MT_Success_' cortical_region '.mat'], 'rt_mt_perform_Info')
% rt_mt_perform_Info(end) = [];

reorder = true; %Note this only changes black avg dots
rTmT_calc = 'median'; % 'mean or 'median'
medPrct = [25 75];

plotInd = false;  %To toggle plotting each target separately on (true) and off (false)

allFig = figure('Position', [400 100 400 800]);
all_lineFS = 16;
all_lineW = 1;
all_axFS = 18;
markerSize = 60;

rt_yLim = [0 1200];
rt_yTicks = [0:300:1200];
mt_yLim = [0 1200];
mt_yTicks = [0:300:1200];

monkSettings = getMonkSettings(monk, cortical_region, rt_mt_perform_Info);
useColor = monkSettings.useColor;
skipNum = monkSettings.skipNum;
dimStart = monkSettings.st_fu_dr_si(1);
dimFull = monkSettings.st_fu_dr_si(2);
dropSession = monkSettings.st_fu_dr_si(3);
singleSession = monkSettings.st_fu_dr_si(4);

if useResT
    numT = 5; %Keep at 5 since using nanmean and including res targets now
else
    numT = 4;
end

% Extract information from saved variable for plotting
fileList = string(nan(length(rt_mt_perform_Info), size(useColor,2)));
dateList = string(nan(length(rt_mt_perform_Info), size(useColor,2)));
plotColors = string(nan(length(rt_mt_perform_Info), size(useColor,2)));
prctColors = zeros(length(rt_mt_perform_Info), size(useColor,2));

meanRT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
meanMT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
stdRT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
stdMT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
medRT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
medMT_all = nan(length(rt_mt_perform_Info), size(useColor,2));
prctRT_all = cell(length(rt_mt_perform_Info), size(useColor,2));
prctMT_all = cell(length(rt_mt_perform_Info), size(useColor,2));
allSuccess = nan(length(rt_mt_perform_Info), size(useColor,2));
allTrials = nan(length(rt_mt_perform_Info), 1);


scatX = cell(6, size(useColor,2)); 
successY = cell(size(scatX));
rtY = cell(size(scatX));
mtY = cell(size(scatX));
ptchX = cell(size(scatX));
rtBars = cell(size(scatX));
mtBars = cell(size(scatX));

for ci = 1:size(useColor,2)
    
    meanRT_targets = nan(length(rt_mt_perform_Info),numT);
    meanMT_targets = nan(length(rt_mt_perform_Info),numT);
    targetSuccess = nan(length(rt_mt_perform_Info),numT);

    for fi = 1:length(rt_mt_perform_Info)
         
        if ~isempty(vertcat(rt_mt_perform_Info(fi).targetRT{:}))
            dateList(fi,ci) = string(rt_mt_perform_Info(fi).Date);
            [thisPath, thisFile, thisExt] = fileparts(rt_mt_perform_Info(fi).File);
            fileList(fi,ci) = string([thisFile thisExt]);
            
            if strcmpi(useColor(fi,ci), "Max") && ~strcmpi(rt_mt_perform_Info(fi).colorFreqs, "NaN")
                colorNDX = rt_mt_perform_Info(fi).colorFreqs == max(rt_mt_perform_Info(fi).colorFreqs);
                plotColors(fi,ci) = rt_mt_perform_Info(fi).allColors(colorNDX);
                prctColors(fi,ci) = max(rt_mt_perform_Info(fi).colorFreqs);
            elseif strcmpi(useColor(fi,ci), "All") || strcmpi(rt_mt_perform_Info(fi).colorFreqs, "NaN")
                colorNDX = true(length(rt_mt_perform_Info(fi).allColors),1);
                plotColors(fi,ci) = "All";
                prctColors(fi,ci) = 1;
            else 
                colorNDX = strcmpi(rt_mt_perform_Info(fi).allColors, useColor(fi,ci));
                plotColors(fi,ci) = useColor(fi,ci);
                prctColors(fi,ci) = rt_mt_perform_Info(fi).colorFreqs(colorNDX);
            end
            
            for targ_i = 1:size(rt_mt_perform_Info(fi).targetRT, 2)
                useRT{targ_i} = vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,targ_i});
                useMT{targ_i} = vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,targ_i});
            end
    
            meanRT_targets(fi,1:size(rt_mt_perform_Info(fi).targetRT, 2)) = cellfun(@(x) nanmean(x), useRT);
            meanMT_targets(fi,1:size(rt_mt_perform_Info(fi).targetMT, 2)) = cellfun(@(x) nanmean(x), useMT);
            targetSuccess(fi,1:size(rt_mt_perform_Info(fi).targetSuccess, 2)) = rt_mt_perform_Info(fi).targetSuccess(colorNDX,:)./rt_mt_perform_Info(fi).targetTrials(colorNDX,:);
            
            meanRT_all(fi,ci) = nanmean(vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,:}));
            stdRT_all(fi,ci) = nanstd(vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,:}));
            meanMT_all(fi,ci) = nanmean(vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,:}));
            stdMT_all(fi,ci) = nanstd(vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,:}));
            allSuccess(fi,ci) = nansum(rt_mt_perform_Info(fi).targetSuccess(colorNDX,:), 'all')/nansum(rt_mt_perform_Info(fi).targetTrials(colorNDX,:), 'all');  
    
            medRT_all(fi,ci) = nanmedian(vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,:}));
            prctRT_all{fi,ci} = [prctile(vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,:}), medPrct(1)), prctile(vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,:}), medPrct(2))];
            medMT_all(fi,ci) = nanmedian(vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,:}));
            prctMT_all{fi,ci} = [prctile(vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,:}), medPrct(1)), prctile(vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,:}), medPrct(2))];
            allSuccess(fi,ci) = nansum(rt_mt_perform_Info(fi).targetSuccess(colorNDX,:), 'all')/nansum(rt_mt_perform_Info(fi).targetTrials(colorNDX,:), 'all');
            allTrials(fi,1) = nansum(rt_mt_perform_Info(fi).targetTrials, 'all');
        
        end
        
        if any(fi == skipNum)   % In case certain days were omitted
            dateList(fi,ci) = NaN;
            fileList(fi,ci) = NaN;
            plotColors(fi,ci) = "NaN";
            prctColors(fi,ci) = NaN;
            meanRT_targets(fi,:) = NaN;
            meanMT_targets(fi,:) = NaN;
            meanRT_all(fi,ci) = NaN;
            meanMT_all(fi,ci) = NaN;
            stdRT_all(fi,ci) = NaN;
            stdMT_all(fi,ci) = NaN;
            medMT_all(fi,:) = NaN;
            medRT_all(fi,ci) = NaN;
            prctMT_all{fi,ci} = NaN;
            prctRT_all{fi,ci} = NaN;
            allSuccess(fi,ci) = NaN;
            allTrials(fi,1) = NaN;
        end  
        
    end
    
    useDays = find(~ismissing(fileList(:,ci)));
    plotDays = [1:length(useDays)]';
    xLim = [1 length(useDays)];
    if reorder
        plotOrder = monkSettings.plotOrder;
    else
        plotOrder = [1:length(useDays)]'; %If you want to plot in a different order from chronological by date, change this 
    end
    
    if useCatchT
        targetStrs = ["Target 0", "Target 2", "Target 4", "Target 6", "Catch Target"];
    elseif useResT
        targetStrs = ["Target 0", "Target 2", "Target 4", "Target 6", "Res Target"];
    else
        targetStrs = ["Target 0", "Target 2", "Target 4", "Target 6"];
    end
    allNum = size(targetStrs,2)+1;
    
    if plotInd
        rtFig = figure('Position', [100 400 550 350]);
        title('Reaction Times')
        xlabel('Session')
        ylabel('Time (msec)')
        ylim([100 1000])
    
        mtFig = figure('Position', [700 400 550 350]);
        title('Movement Times')
        xlabel('Session')
        ylabel('Time (msec)')
        ylim([100 600])
    
        successFig = figure('Position', [1300 400 550 350]);
        title('Success Percentages')
        xlabel('Session')
        ylabel('Percent Correct')
        ylim([0 100])
    
        for ti = 1:size(targetStrs,2)
    
            figure(rtFig)
            hold on
            rtScat(ti) = scatter(plotDays, 1000*meanRT_targets(useDays(plotOrder),ti), markerSize/2, 'filled');
            rtScat(ti).DisplayName = targetStrs(ti);
            hold off
    
            figure(mtFig)
            hold on
            mtScat(ti) = scatter(plotDays, 1000*meanMT_targets(useDays(plotOrder),ti), markerSize/2, 'filled');
            mtScat(ti).DisplayName = targetStrs(ti);
            hold off
    
            figure(successFig)
            hold on
            successScat(ti) = scatter(plotDays, 100*targetSuccess(useDays(plotOrder),ti), markerSize/2, 'filled');
            successScat(ti).DisplayName = targetStrs(ti);
            hold off
    
        end
        
        figure(rtFig)
        hold on
        rtScat(allNum) = scatter(plotDays, 1000*meanRT_all(useDays(plotOrder)), markerSize, 'k', 'filled');
        rtScat(allNum).DisplayName = 'All Targets';
        if ~isnan(dropSession)
            rtScat(allNum+1) = xline(dropSession, '--k', 'Single Electrode', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        end
        rtScat(allNum+2) = xline(dimStart, '--g', 'Dim Start', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        rtScat(allNum+3) = xline(dimFull, '--r', 'ICMS Only', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        hold off    
        legend(rtScat(1:allNum), 'location', 'northwest')    
    
        figure(mtFig)
        hold on
        mtScat(allNum) = scatter(plotDays, 1000*meanMT_all(useDays(plotOrder)), markerSize, 'k', 'filled');
        mtScat(allNum).DisplayName = 'All Targets';
        if ~isnan(dropSession)
            mtScat(allNum+1) = xline(dropSession, '--k', 'Single Electrode', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        end
        mtScat(allNum+2) = xline(dimStart, '--g', 'Dim Start', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        mtScat(allNum+3) = xline(dimFull, '--r', 'ICMS Only', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        hold off        
        legend(mtScat(1:allNum), 'location', 'northwest' )
    
        figure(successFig)
        hold on
        successScat(allNum) = scatter(plotDays, 100*allSuccess(useDays(plotOrder)), markerSize, 'k', 'filled');
        successScat(allNum).DisplayName = 'All Targets';
        if ~isnan(dropSession)
            successScat(allNum+1) = xline(dropSession, '--k', 'Single Electrode', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        end
        successScat(allNum+2) = xline(dimStart, '--g', 'Dim Start', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        successScat(allNum+3) = xline(dimFull, '--r', 'ICMS Only', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        hold off        
        legend(successScat(1:allNum), 'location', 'southwest')
        
    end
    
    
    %%% Combination Plot %%%
    figure(allFig)
    
    dimNoneNDX = (1:length(plotOrder))' < dimStart;
    dropNDX = (1:length(dimNoneNDX))' >= dropSession & (1:length(dimNoneNDX))' <= singleSession;
    singleNDX = (1:length(dimNoneNDX))' >= singleSession;
    dimPaleNDX = strcmpi(plotColors(useDays(plotOrder),ci), monkSettings.dimStr);% | strcmpi(plotColors(useDays(plotOrder),ci), "dim5d");
    dimFullNDX = strcmpi(plotColors(useDays(plotOrder),ci), "dimFull") & ~dropNDX & ~singleNDX;
    dimSomeNDX = ~dimNoneNDX & ~dimFullNDX & ~dimPaleNDX & ~dropNDX & ~singleNDX;
    
    levelNDX = [dimNoneNDX,  dimSomeNDX,  dimPaleNDX,  dimFullNDX,  dropNDX,  singleNDX];
    if strcmpi(monk(1), 'F') && strcmpi(cortical_region, 'S1')
        levelNDX(35,:) = [false, false, false, true, false, false];
    end
    if strcmpi(monk(1), 'Q') && strcmpi(cortical_region, 'Pmv')
        levelNDX(18,:) = [false, false, false, true, false, false];
    end

    for level = 1:size(levelNDX,2)
        
        thisLevel = levelNDX(:,level);
        scatX{level,ci} = plotDays(thisLevel);
        successY{level,ci} = 100*allSuccess(useDays(plotOrder(thisLevel)),ci);
        if sum(thisLevel) > 0 
            switch rTmT_calc
                case 'mean'
                    rtY{level,ci} = 1000*meanRT_all(useDays(plotOrder(thisLevel)),ci);
                    mtY{level,ci} = 1000*meanMT_all(useDays(plotOrder(thisLevel)),ci);
            
                    ptchX{level,ci} = [plotDays(thisLevel), flip(plotDays(thisLevel))]';
                    rtBars{level,ci} = 1000*[(meanRT_all(useDays(plotOrder(thisLevel)),ci)-stdRT_all(useDays(plotOrder(thisLevel)),ci)); flip((meanRT_all(useDays(plotOrder(thisLevel)),ci)+stdRT_all(useDays(plotOrder(thisLevel)),ci)))];
                    mtBars{level,ci} = 1000*[(meanMT_all(useDays(plotOrder(thisLevel)),ci)-stdMT_all(useDays(plotOrder(thisLevel)),ci)); flip((meanMT_all(useDays(plotOrder(thisLevel)),ci)+stdMT_all(useDays(plotOrder(thisLevel)),ci)))];
                case 'median'
                    rtY{level,ci} = 1000*medRT_all(useDays(plotOrder(thisLevel)),ci);
                    mtY{level,ci} = 1000*medMT_all(useDays(plotOrder(thisLevel)),ci);
            
                    ptchX{level,ci} = [plotDays(thisLevel); flip(plotDays(thisLevel))]; 
                    thisPrctRT = vertcat(prctRT_all{useDays(plotOrder(thisLevel)),ci}); 
                    thisPrctMT = vertcat(prctMT_all{useDays(plotOrder(thisLevel)),ci});
                    rtBars{level,ci} = 1000*[thisPrctRT(:,1); flip(thisPrctRT(:,2))];
                    mtBars{level,ci} = 1000*[thisPrctMT(:,1); flip(thisPrctMT(:,2))];
            end
        end
    end
end

 
levelColors = {'k', [0 0.6 0], [0.4 1 0.4], 'r', 'b', [0.6 0.6 0.6]};
% levelColors = {'k', 'k', 'k'};
pltLevels = [1 1 1 1 1 1; 1 1 2 2 1 1];   % Here, plt Level refers to the difficulty when 2 levels are present during a single day (i.e dim5c and Full)
markerSize = 120;
for di = [2 1]
    for level = [1 2 4 3 5 6]%1:size(levelNDX,2)
    
        pltNDX = pltLevels(di,level);
        if pltNDX == 2
            isUnique = ~ismember(scatX{level,2}, scatX{level,1});  
            if ~any(isUnique)
                continue
            end
            isUnique(max(find(isUnique))+1) = true;
        elseif di == 2
            continue
        else
            isUnique = true(size(scatX{level,1}));
        end
        if isempty(scatX{level,di})
            continue
        end
    
        subplot(3,1,1)
        ax = gca; ax.LineWidth = 2;
        hold on
        successScat(level, di) = scatter(scatX{level,pltNDX}, successY{level,pltNDX}, markerSize, 'filled', 'MarkerFaceColor', levelColors{level});
        successScat(level, di).DisplayName = 'All Targets';
        if ~isnan(dropSession)
            successScat(level+1, di) = xline(dropSession, '--b', 'LineWidth', 2);%, 'Drop', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        if ~isnan(singleSession)
            successScat(level+2, di) = xline(singleSession, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);%, 'Single', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        successScat(level+3, di) = xline(dimStart, '--', 'Color', [0 0.7 0.1], 'LineWidth', 2);%, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        successScat(level+4, di) = xline(dimFull, '--r', 'LineWidth', 2);%, 'ICMS', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        hold off        
        ylabel('% Correct')
        yticks(0:20:100)
        yticklabels('')
        ylabel('')
        ylim([0 100])
        xlim(xLim) 
        set(gca, 'FontSize', all_axFS) 
        if useTitle
            title([monk ', ' cortical_region], 'FontSize', 14)
        end
        
        if strcmpi(cortical_region, 'AIP') && strcmpi(monk(1), 'F')
    %         xline(27, '--', 'Rem T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(36, '--', 'Drop', 'Color', 'b', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(37, '--', 'Single', 'Color', 'k', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(40, '--', 'Add T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
            xline(26, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(31, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(36, '--', '', 'Color', 'b', 'LineWidth', 2);
            xline(37, '--', '', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            xline(40, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'Q')
            xline(37, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'F')
            xline(20, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(32, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(39, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        
        
        subplot(3,1,2)
        ax = gca; ax.LineWidth = 2;
        hasSuccess = ~isnan(rtBars{level,pltNDX});
        hold on
        rtPatch(level, di) = patch(ptchX{level,pltNDX}(hasSuccess&[isUnique; flip(isUnique)]), rtBars{level,pltNDX}(hasSuccess&[isUnique; flip(isUnique)]), levelColors{level}, 'edgecolor', 'none', 'FaceAlpha', 0.2);
        if ~isnan(dropSession)
            rtScat(level+1, di) = xline(dropSession, '--b', 'LineWidth', 2);%, 'Drop', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        if ~isnan(singleSession)
            rtScat(level+2, di) = xline(singleSession, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);%, 'Single', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        rtScat(level+3, di) = xline(dimStart, '--', 'Color', [0 0.7 0.1], 'LineWidth', 2);%, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        rtScat(level+4, di) = xline(dimFull, '--r', 'LineWidth', 2);%, 'ICMS', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        rtScat(level, di) = scatter(scatX{level,pltNDX}, rtY{level,pltNDX}, markerSize, 'filled', 'MarkerFaceColor', levelColors{level});
        rtScat(level, di).DisplayName = 'All Targets';
        hold off   
        ylabel('RT (msec)')
        yticks(rt_yTicks)
        yticklabels('')
        ylabel('')
        ylim(rt_yLim)
        xlim(xLim)
        set(gca, 'FontSize', all_axFS)
        
        if strcmpi(cortical_region, 'AIP') && strcmpi(monk(1), 'F')
    %         xline(27, '--', 'Rem T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(36, '--', 'Drop', 'Color', 'b', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(37, '--', 'Single', 'Color', 'k', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(40, '--', 'Add T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
            xline(26, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(31, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(36, '--', '', 'Color', 'b', 'LineWidth', 2);
            xline(37, '--', '', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            xline(40, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'Q')
            xline(37, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'F')
            xline(20, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(32, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(39, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        
        subplot(3,1,3)
        ax = gca; ax.LineWidth = 2;
        hold on
        mtPatch(level, di) = patch(ptchX{level,pltNDX}(hasSuccess&[isUnique; flip(isUnique)]), mtBars{level,pltNDX}(hasSuccess&[isUnique; flip(isUnique)]), levelColors{level}, 'edgecolor', 'none', 'FaceAlpha', 0.2);
        if ~isnan(dropSession)
            mtScat(level+1, di) = xline(dropSession, '--b', 'LineWidth', 2);%, 'Drop', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        if ~isnan(singleSession)
            mtScat(level+2, di) = xline(singleSession, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);%, 'Single', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        end
        mtScat(level+3, di) = xline(dimStart, '--', 'Color', [0 0.7 0.1], 'LineWidth', 2);%, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        mtScat(level+4, di) = xline(dimFull, '--r', 'LineWidth', 2);%, 'ICMS', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
        mtScat(level, di) = scatter(scatX{level,pltNDX}, mtY{level,pltNDX}, markerSize, 'filled', 'MarkerFaceColor', levelColors{level});
        mtScat(level, di).DisplayName = 'All Targets';
        hold off    
        xlabel('Session')
        ylabel('MT (msec)')
        yticks(mt_yTicks)
        yticklabels('')
        ylabel('')
        ylim([mt_yLim])
        xlim(xLim)
        set(gca, 'FontSize', all_axFS)
        % legend(mtScat(allNum), 'location', 'northeast')
        
        if strcmpi(cortical_region, 'AIP') && strcmpi(monk(1), 'F')
    %         xline(27, '--', 'Rem T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(36, '--', 'Drop', 'Color', 'b', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(37, '--', 'Single', 'Color', 'k', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
    %         xline(40, '--', 'Add T4', 'Color', [253, 142, 51]/255, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', 'FontSize', all_lineFS, 'LineWidth', all_lineW);
            xline(26, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(31, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(36, '--', '', 'Color', 'b', 'LineWidth', 2);
            xline(37, '--', '', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            xline(40, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'Q')
            xline(37, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
        if strcmpi(cortical_region, 'S1') && strcmpi(monk(1), 'F')
            xline(20, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(32, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
            xline(39, '--', '', 'Color', [253, 142, 51]/255, 'LineWidth', 2);
        end
    end
end

useFiles = fileList(useDays,1); useColors = plotColors(useDays,:); usePrct = prctColors(useDays,:);
theseFiles = [useDays(plotOrder), plotDays, useFiles(plotOrder)]
trialCounts = round([nanmean(allTrials), nanstd(allTrials)])
% trialCounts = round([nanmean(allTrials), min(allTrials), max(allTrials)])

rtInfo = rt_mt_perform_Info;
rtSessions = [21:24, 28];
dimColor = "dimFull";
saveRTs(monk, cortical_region, rtSessions, dimColor, rtInfo)
saveas(gcf, ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\' 'SP_RT_MT_OverTime.png'])
saveas(gcf, [savePath  monk '_' cortical_region '_' 'SP_RT_MT_OverTime.png'])
print([savePath  monk '_' cortical_region '_' 'SP_RT_MT_OverTime.tiff'], '-dtiff', '-r300')
% [plotDays(thisLevel)', fileList, allSuccess, meanRT_all, meanMT_all]
% [plotDays(thisLevel)', fileList]
% [useDays(plotOrder), plotDays, useColors(plotOrder,:), usePrct(plotOrder,:)]





%% Extra functions

function out = getTrialID(trialStart, markers)
    
    if ~isempty(trialStart)
        for ti = 1:length(trialStart)
            out(ti,:) = markers(trialStart(ti)+1:trialStart(ti)+4);
        end
    else
        out = NaN;
    end
    
end

function monkSettings = getMonkSettings(monk, cortical_region, rt_mt_perform_Info)
    
    useColor = string(nan(length(rt_mt_perform_Info),2));
    useColor(:,:) = "Max";
    switch monk(1)
        case 'Q'
            dimStr = "dim5d";
            switch cortical_region
                case 'S1'
                    skipNum = [58];
                    if length(rt_mt_perform_Info) > skipNum
                        numDays = length(rt_mt_perform_Info) - length(skipNum);
                    else
                        numDays = length(rt_mt_perform_Info);
                    end
                    plotOrder = [1:numDays]';

                    st_fu_dr_si = [6, 28, 50, 51];
                    useColor(28:end,:) = "dimFull";
                    useColor(22:33,1) = "dim5d";                        
                case 'PMv'
                    skipNum = [34];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:10, 12, 11, 13:numDays]';

                    st_fu_dr_si = [5, 14, 27, 29];
                    useColor(14:end,:) = "dimFull";
                    useColor([11:14],1) = "dim5d"; 
                case 'AIP'
                    skipNum = [6, 8, 45];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:numDays]';

                    st_fu_dr_si = [12, 21, 28, 30];
                    useColor(23:end,:) = "dimFull";
                    useColor(23,1) = "dim5d"; 
                case 'PMd'
                    skipNum = [10, 11, 12, 32];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:24, 28, 25:27, 29:numDays]';
                    % plotOrder = [1:numDays]';

                    st_fu_dr_si = [10, 20, NaN, NaN]; 
                    useColor([24:end],:) = "dimFull";
                    useColor(18:26,1) = "dim5d";
                case 'dPPC'
                    skipNum = []; 
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:length(rt_mt_perform_Info)]';

                    st_fu_dr_si = [10, 16, NaN, NaN]; 
                    useColor(16:end,:) = "dimFull";
                    useColor(16:19,1) = "dim5d";
            end
        case 'F'
            dimStr = "dim5e";
            switch cortical_region
                case 'S1'
                    skipNum = [4, 40]; 
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    % plotOrder = [1:11, 18, 15:17, 12:14, 19:numDays]';  
                    plotOrder = [1:numDays]'; 

                    st_fu_dr_si = [4, 12, 43, 45];
                    useColor(13:15,1) = "dim5d";
                    useColor([16,18,19,20],1) = "dim5e";
                    
                    useColor([18],2) = "dim5e";
                    useColor([13:15, 19,20],2) = "dimFull";
                    % useColor([13:14],1) = "dim5d";
                    % useColor([16:19],:) = "dim5e";
                    x=2;
                case 'PMv'
                    skipNum = [1, 2, 3];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:numDays]';

                    st_fu_dr_si = [4, 11, 19, 21];
                    useColor(14:end,:) = "dimFull";
                    useColor(13:15,1) = "dim5d";
%                     useColor(14,1) = "dim5e";
                case 'AIP'
                    skipNum = [58];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
%                     plotOrder = [1:numDays]';
                    plotOrder = [1:16, 20, 17:19, 21:23, 25:26, 24, 27:numDays]; 

                    st_fu_dr_si = [10, 17, 46, 47];
                    useColor(17:end,:) = "dimFull";
                    useColor([16,17,20],1) = "dim5e";
                case 'PMd'
                    skipNum = [1:5, 27,28];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:numDays]';

                    st_fu_dr_si = [10, 18, NaN, NaN];
                    useColor(23:end,:) = "dimFull";
                    useColor(21:25,1) = "dim5e";
                    useColor(28,:) = "dim5e";
                case 'dPPC'
                    skipNum = [29];
                    numDays = length(rt_mt_perform_Info) - length(skipNum);
                    plotOrder = [1:24, 28, 25:27]';

                    st_fu_dr_si = [3, 21, NaN, NaN]; 
                    useColor(21:end,:) = "dimFull";
                    useColor(19:24,1) = "dim5d";
            end
    end
   
    monkSettings.useColor = useColor;
    monkSettings.skipNum = skipNum;
    monkSettings.st_fu_dr_si = st_fu_dr_si;
    monkSettings.plotOrder = plotOrder;
    monkSettings.dimStr = dimStr;

end


function saveRTs(monk, cortical_region, rtSessions, dimColor, rtInfo)
    
    matFile = 'F:\Publications\InjectInformationPMCPPC\Results\rtSuccessComp.mat';
    if exist(matFile, 'file')
        load('F:\Publications\InjectInformationPMCPPC\Results\rtSuccessComp.mat', 'rtSuccessStruct')
        rtSuccessStruct(end+1).Monk = monk;

    else
        rtSuccessStruct = struct();
    end
    % rtSessions = [11 12 13 14];
    % dimColor = "dim5d";
    theseRTs = [];
    theseTrials = [];
    theseSuccess = [];
    for si = rtSessions
    
        clrNDX = strcmpi(rtInfo(si).allColors, dimColor);% | strcmpi(rtInfo(si).allColors, "dim5e");
        theseRTs = [theseRTs; vertcat(rtInfo(si).targetRT{clrNDX,:})];
        theseTrials = [theseTrials; vertcat(rtInfo(si).targetTrials(clrNDX,:))];
        theseSuccess = [theseSuccess; vertcat(rtInfo(si).targetSuccess(clrNDX,:))];
    
    end
    rtSuccessStruct(end).Monk = monk;
    rtSuccessStruct(end).Area = cortical_region;
    rtSuccessStruct(end).rts = theseRTs;
    rtSuccessStruct(end).trialNums = theseTrials;
    rtSuccessStruct(end).successNums = theseSuccess;
    save('F:\Publications\InjectInformationPMCPPC\Results\rtSuccessComp.mat', 'rtSuccessStruct')
end
