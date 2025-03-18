%% Add paths to necessary toolboxes

clear; close all; 
addpath(genpath('D:\Projects\COT_ICMS\distantMod\'))
addpath(genpath('D:\MATLAB\FML_matlab\'))
addpath(genpath('D:\MATLAB\fileHandling_BR\'))
addpath(genpath('D:\MATLAB\Statistics\'))

%% Call Files to analyze

monk ='Qulio';
cortical_region = 'AIP';

useResT = false;
useCatchT = false;
errorBars = 'median'; % 'mean or 'median'
medPrct = [25 75];

ResT = 3;
CatchT = 3;

path_to_data = ['D:\RawData\COT_ICMS\' monk '\' cortical_region '\'];
savePath = ['D:\Publications\InjectInformationPMCPPC\Results\'];
all_nevFiles = [dir([path_to_data 'Ripple\' monk(1) '_*ICMS_*.nev']); dir([path_to_data 'Ripple\paramSweeps\' monk(1) '_*ICMS_*.nev'])];
all_ns2Files = [dir([path_to_data 'Ripple\' monk(1) '_*ICMS_*.ns2']); dir([path_to_data 'Ripple\paramSweeps\' monk(1) '_*ICMS_*.ns2'])];
all_h5Files = [dir([path_to_data 'UTCS\' monk(1) '_*ICMS_*.h5']); dir([path_to_data 'UTCS\paramSweeps\' monk(1) '_*ICMS_*.h5'])];
nev_fileList = {all_nevFiles(:).name}'; nev_folderList = {all_nevFiles(:).folder}';
ns2_fileList = {all_ns2Files(:).name}'; ns2_folderList = {all_ns2Files(:).folder}';
h5_fileList = {all_h5Files(:).name}'; h5_folderList = {all_h5Files(:).folder}';

% Check for existing .mat file with previous information and list of used files
% If there is no .mat file a new one will be created

if exist([path_to_data monk '_all_rTmTsuccess_' cortical_region '.mat'], 'file') ~= 0 
    load([path_to_data monk '_all_rTmTsuccess_' cortical_region '.mat'], 'rt_mt_perform_Info')           
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
nev_fileList = nev_fileList(newNDX); nev_folderList = nev_folderList(newNDX);
ns2_fileList = ns2_fileList(newNDX); ns2_folderList = ns2_folderList(newNDX);
h5_fileList = h5_fileList(newNDX); h5_folderList = h5_folderList(newNDX);
[ns2_fileList, nev_fileList, h5_fileList]
correctEventCodes = false(size(ns2_fileList));
if strcmpi(monk(1), 'Q')
    correctEventCodes(cellfun(@(x) contains(x, '20210413')|contains(x, '20210419'), ns2_fileList)) = true;
    correctEventCodes(cellfun(@(x) contains(x, '20211015')|contains(x, '20211018'), ns2_fileList)) = true;
end
numUsedFiles = sum(~cellfun(@(x) isempty(x), {rt_mt_perform_Info(:).File}'));
%%

startTime = tic;
for fi = 1:length(ns2_fileList)
    
    disp(['File ' num2str(fi) ' of ' num2str(sum(newNDX)) ' new files'])
    rt_mt_NDX = fi+numUsedFiles;  
    joystick_data_filename = [ns2_folderList{fi} '\' ns2_fileList{fi}];
    marker_data_filename = [nev_folderList{fi} '\' nev_fileList{fi}];
    h5_filename = [h5_folderList{fi} '\' h5_fileList{fi}];
    
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
        h5_filename = [h5_folderList{fi} '\' h5_fileList{fi}];
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
        for i=1:trials_count
            ec{i} = h5read(h5_filename, ['/trial_' pad(num2str(i), 4, 'left', '0'), '/event_codes']);
        end
        [clrs, changes]= COTReportTargetColorChange(h5_filename);
        clrs = translateClrs(clrs);
        theseClrs = unique(clrs);
        existClrsAtt = true;
    catch
        existClrsAtt = false;
        warning('No Target Color Attribute')
        theseClrs = "No Color Attribute in .h5";
    end

    if length(js_sp) > 1
        allFields = {'baseline_mean'; 'baseline_stddev'; 'onset_time'; 'reaction_duration'; 'movement_duration'; 'mv_onset'};
        for frag = 1:length(js_sp)
             
            % extract trials and align on first marker
            [js_T, js_MT, js_patsIdx, js_patsT] = TimeseriesFindExtractAlign(js_sp(frag), js_markers_ts, js_markers_sv, successPatterns);
            
            % calculate movement onset
            js_onsets{frag,1} = FMLFindOnset(js_T, js_MT, 'BaselineFromIndex', 2);
            
            % to modify FMLFindOnset options:
            % onset_options = FMLFindOnset; % default options
            % onset_options.baselinefromindex = 2;  % index of marker that anchors the baseline calculation
            % onset_options.baselineduration = [-0.5, 0]; % from 500ms before baseline marker time, to 0ms before the baseline marker
            % onset_options.onsetfromindex = 3;    % index of the marker from which time we seach for onset
            % onset_options.onsetstddevfactor = 3; % onset is when over this factor times stddev of baseline from baseline...
            % onset_options.onsetduration = 0.03; % ... from this amount of time
            % onset_options.reachtoindex = 5; % to calculate movement duration: from onset to this marker
            % js_onsets = FMLFindOnset(js_T, js_MT, onset_options);
    
            for targ_i = 1:length(js_T)
                fragID = cellfun(@(x) isempty(x.Data), js_T{targ_i});
                js_onsets{frag}{targ_i}(~fragID) = [];
            
                missingFields = cellfun(@(x) ~isfield(js_onsets{frag}{targ_i}, x), allFields);
                for thisField = 1:length(missingFields)
                    if missingFields(thisField)
                        js_onsets{frag}{targ_i} = setfield(js_onsets{frag}{targ_i}, allFields{thisField}, []);
                    end
                end
    
            end 
        end
    
        for frag = 2:length(js_onsets)
            for targ_i = 1:length(js_onsets{1})
                js_onsets{1}{targ_i} = [js_onsets{1}{targ_i}, js_onsets{frag}{targ_i}];
            end
        end
        js_onsets = js_onsets{1};
    else
        [js_T, js_MT, js_patsIdx, js_patsT] = TimeseriesFindExtractAlign(js_sp, js_markers_ts, js_markers_sv, successPatterns);
        js_onsets = FMLFindOnset(js_T, js_MT, 'BaselineFromIndex', 2);
    end

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

%Add new files into existing list chronilogically by date

[~, sortNDX] = sort(allDates);
rt_mt_perform_Info = rt_mt_perform_Info(sortNDX);

save([path_to_data monk '_all_rTmTsuccess_' cortical_region '.mat'], 'rt_mt_perform_Info') %Save information


%% Performance Calculations Data

monks =  ["Qulio", "Felix"];
corticalAreas = ["S1", "PMv", "AIP"];
rTmT_calc = 'median';
medPrct = [25 75];

for mi = 1:length(monks)
    for ci = 1:length(corticalAreas)

        monk = char(monks(mi));
        cortical_region = char(corticalAreas(ci));
        path_to_data = ['D:\RawData\COT_ICMS\' monk '\' cortical_region '\'];
    
        load([path_to_data monk '_all_rTmTsuccess_' cortical_region '.mat'], 'rt_mt_perform_Info')

        switch monk(1)
            case 'Q'
                switch cortical_region
                    case 'S1'
                        skipNum = [58];
                    case 'PMv'
                        skipNum = [34, 35];           % PMv
                    case 'AIP'
                        skipNum = [6, 8, 45, 46];
                end
            case'F'
                switch cortical_region
                    case 'S1'
                        skipNum = [4, 51]; 
                    case 'PMv'
                        skipNum = [1, 2, 3, 29];
                    case 'AIP'
                        skipNum = [58];
                end
        end 
        rt_mt_perform_Info(skipNum) = [];
    
        useColor = string(nan(length(rt_mt_perform_Info),1));
        useColor(:) = "Max";
        switch monk(1)
            case 'Q'
                switch cortical_region
                    case 'S1'
                        strt_end = [5, 51];                         
                    case 'PMv'
                        strt_end = [5, 29];        
                    case 'AIP'
                        strt_end = [12, 30];                        
                end
            case'F'
                switch cortical_region
                    case 'S1'
                        strt_end = [4, 45];
                    case 'PMv'
                        strt_end = [4, 21];  
                    case 'AIP'
                        strt_end = [10, 47];
                        useColor(17) = "dimFull";
                        useColor(20:end) = "dimFull";  
                end
        end
    
        for fi = 1:length(rt_mt_perform_Info)
             
            if ~isempty(vertcat(rt_mt_perform_Info(fi).targetRT{:}))
                dateList(fi, 1) = string(rt_mt_perform_Info(fi).Date);
                [thisPath, thisFile, thisExt] = fileparts(rt_mt_perform_Info(fi).File);
                fileList(fi, 1) = string([thisFile thisExt]);
                
                if strcmpi(useColor(fi), "Max") && ~strcmpi(rt_mt_perform_Info(fi).colorFreqs, "NaN")
                    colorNDX = rt_mt_perform_Info(fi).colorFreqs == max(rt_mt_perform_Info(fi).colorFreqs);
                    plotColors(fi) = rt_mt_perform_Info(fi).allColors(colorNDX);
                    prctColors(fi) = max(rt_mt_perform_Info(fi).colorFreqs);
                elseif strcmpi(useColor(fi), "All") || strcmpi(rt_mt_perform_Info(fi).colorFreqs, "NaN")
                    colorNDX = true(length(rt_mt_perform_Info(fi).allColors),1);
                    plotColors(fi) = "All";
                    prctColors(fi) = 1;
                else 
                    colorNDX = strcmpi(rt_mt_perform_Info(fi).allColors, useColor(fi));
                    plotColors(fi) = useColor(fi);
                    prctColors(fi) = rt_mt_perform_Info(fi).colorFreqs(colorNDX);
                end
                
                for targ_i = 1:size(rt_mt_perform_Info(fi).targetRT, 2)
                    useRT{fi,targ_i} = vertcat(rt_mt_perform_Info(fi).targetRT{colorNDX,targ_i});
                    useMT{fi,targ_i} = vertcat(rt_mt_perform_Info(fi).targetMT{colorNDX,targ_i});
                end
        
                meanRT_targets(fi,1:size(rt_mt_perform_Info(fi).targetRT, 2)) = cellfun(@(x) nanmean(x), useRT);
                meanMT_targets(fi,1:size(rt_mt_perform_Info(fi).targetMT, 2)) = cellfun(@(x) nanmean(x), useMT);
                targetSuccess(fi,1:size(rt_mt_perform_Info(fi).targetSuccess, 2)) = rt_mt_perform_Info(fi).targetSuccess(colorNDX,:);
                targetNum(fi,1:size(rt_mt_perform_Info(fi).targetSuccess, 2)) = rt_mt_perform_Info(fi).targetTrials(colorNDX,:);
               
            
            end
        end

%         figure()
%         plot(100*allSuccess, 'o', 'MarkerFaceColor', 'k')
%         ylim([0 100])
%         xline(strt_end(1))
%         xline(strt_end(2), '--')

        preTrain_success = allSuccess(1:strt_end(1));
        postTrain_success = allSuccess(strt_end(2):end);
        [~, preTrainNDX] = sort(preTrain_success, 'descend');
        [~, postTrainNDX] = sort(postTrain_success, 'descend');

        [hSuccess(mi, ci), pSuccess(mi, ci)] = ttest(preTrain_success(preTrainNDX(1:4)), postTrain_success(postTrainNDX(1:4)));

        switch rTmT_calc
            case 'mean'
                preTrain_RT = meanRT_all(1:strt_end(1));
                postTrain_RT = meanRT_all(strt_end(2):end);
                preTrain_MT = meanMT_all(1:strt_end(1));
                postTrain_MT = meanMT_all(strt_end(2):end);

                preRT(mi,ci) = 1000*mean(preTrain_RT(preTrainNDX(1:4)));
                postRT(mi,ci) = 1000*mean(postTrain_RT(postTrainNDX(1:4)));
                preMT(mi,ci) = 1000*mean(preTrain_MT(preTrainNDX(1:4)));
                postMT(mi,ci) = 1000*mean(postTrain_MT(postTrainNDX(1:4))); 

                [hRT(mi, ci), pRT(mi, ci)] = ttest(preTrain_RT(preTrainNDX(1:4)), postTrain_RT(postTrainNDX(1:4)));
                [hMT(mi, ci), pMT(mi, ci)] = ttest(preTrain_MT(preTrainNDX(1:4)), postTrain_MT(postTrainNDX(1:4)));
            case 'median'
                preTrain_RT = medRT_all(1:strt_end(1));
                postTrain_RT = medRT_all(strt_end(2):end);
                preTrain_MT = medMT_all(1:strt_end(1));
                postTrain_MT = medMT_all(strt_end(2):end);

                preRT(mi,ci) = 1000*median(preTrain_RT(preTrainNDX(1:4)));
                postRT(mi,ci) = 1000*median(postTrain_RT(postTrainNDX(1:4)));
                preMT(mi,ci) = 1000*median(preTrain_MT(preTrainNDX(1:4)));
                postMT(mi,ci) = 1000*median(postTrain_MT(postTrainNDX(1:4)));  
                
                [pRT(mi, ci), hRT(mi, ci)] = signrank(preTrain_RT(preTrainNDX(1:4)), postTrain_RT(postTrainNDX(1:4)));
                [pMT(mi, ci), hMT(mi, ci)] = signrank(preTrain_MT(preTrainNDX(1:4)), postTrain_MT(postTrainNDX(1:4)));

        end       
        
    end
end

preTrain_success
postTrain_success
pSuccess

preRT
postRT
pRT

preMT
postMT
pMT

%%

alpha = 0.05;

% [p, h, stats] = signrank([preSuccess(1,:), preSuccess(2,:)]', [postSuccess(1,:), postSuccess(2,:)]', 'Alpha', alpha)
[h, p, stats] = ttest([preSuccess(1,:), preSuccess(2,:)]', [postSuccess(1,:), postSuccess(2,:)]', 'Alpha', alpha)
figure()
hold on
scatter([preSuccess(1,:), preSuccess(2,:)]', [postSuccess(1,:), postSuccess(2,:)]', 75, 'k')
plot(0:100, 0:100, '--k')
axis square
axis([0 100 0 100])
xticks(0:20:100)
yticks(0:20:100)
xlabel('Post Success')
ylabel('Pre Success')
title(['p = ' num2str(p)])
hold off


% [p, h, stats] = signrank([preRT(1,:), preRT(2,:)]', [postRT(1,:), postRT(2,:)]', 'Alpha', alpha)
[h, p, stats] = ttest([preRT(1,:), preRT(2,:)]', [postRT(1,:), postRT(2,:)]', 'Alpha', alpha)
figure()
hold on
scatter([preRT(1,:), preRT(2,:)]', [postRT(1,:), postRT(2,:)]', 75, 'k')
plot(100:500, 100:500, '--k')
axis square
axis([100 500 100 500])
xticks(100:50:500)
yticks(100:50:500)
xlabel('Post RT')
ylabel('Pre RT')
title(['p = ' num2str(p)])
hold off

% [p, h, stats] = signrank([preMT(1,:), preMT(2,:)]', [postMT(1,:), postMT(2,:)]', 'Alpha', alpha)
[h, p, stats] = ttest([preMT(1,:), preMT(2,:)]', [postMT(1,:), postMT(2,:)]', 'Alpha', alpha)
figure()
hold on
scatter([preMT(1,:), preMT(2,:)]', [postMT(1,:), postMT(2,:)]', 75, 'k')
plot(100:500, 100:500, '--k')
axis square
axis([100 500 100 500])
xticks(100:50:500)
yticks(100:50:500)
xlabel('Post MT')
ylabel('Pre MT')
title(['p = ' num2str(p)])
hold off



%% Extra Functions
function out = getTrialID(trialStart, markers)
    
    if ~isempty(trialStart)
        for ti = 1:length(trialStart)
            out(ti,:) = markers(trialStart(ti)+1:trialStart(ti)+4);
        end
    else
        out = NaN;
    end
    
end

