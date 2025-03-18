clc;
% close all;
loadFiles = logical(1);
if loadFiles
     clear;
     loadFiles = true;
end

%% Adding paths to FML Toolboxes
addpath(genpath('F:\MATLAB\FML_matlab\'))
addpath('F:\MATLAB\fileHandling_BR\')

%% Choose File and sorting parameters

monk = 'Qulio';
cortical_region = 'PMd';
date = '20220209';
dimInfo = 'dim5d_Full_3T';  %Filename after COT_ICMS_
matrixDim = "dimFull"; %Which dim level to show in confusion matrix.
addPath = ''; %'Extra\'
% F_20220622_COT_ICMS_dimFull
% F_20220624_COT_ICMS_dim5e_Full
% F_20220727_COT_ICMS_dimFull_training
% F_20220830_COT_ICMS_dimFull
% F_20220901_COT_ICMS_dimFull_catch


%"dim3", ,"dim4", dim5a", "dim5b", "dim5c", "dim5d", dim5e", "dimFull", or "all" are possible options 
% useBlock = "Full";
useBlock = "exe";
correctEventCodes = false;
usePriorSuccess = true; 
trialsUse = [1:10000]; %Choose trial range or use NaN for all trials

data_path = ['F:\RawData\COT_ICMS\' monk '\' cortical_region '\'];
% data_path = ['E:\QulioBackup\COT_ICMS\S1\'];
addpath(genpath(data_path));
save_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\ConfusionMats\'];  
titleStr = char(matrixDim);
                     
                   
%% Extract trial information from .h5 and .plx file
    % note that the following scripts were designed for UTCS_version 1.0.6
    % changes in .h5 structure may apply
    
numTargets = 8; 
if strcmpi(monk, 'Qulio')
    numTargets = 8;
elseif strcmpi(monk, 'Felix') && str2num(date) >= 20220104
    numTargets = 4;
elseif strcmpi(monk, 'Indira')
    numTargets = 4;
end

if loadFiles

    plx_filename = [data_path 'Plexon\' addPath monk(1) '_' date '_COT_ICMS_' dimInfo '.plx'];
    nev_filename = [data_path 'Ripple\' addPath monk(1) '_' date '_COT_ICMS_' dimInfo '.nev'];
    ns2_filename = [data_path 'Ripple\' addPath monk(1) '_' date '_COT_ICMS_' dimInfo '.ns2'];
    h5_filename = [data_path 'UTCS\' addPath monk(1) '_' date '_COT_ICMS_' dimInfo '.h5'];

%     plx_filename = [data_path 'Plexon\' monk(1) '_' date '_ICMS_' dimInfo '.plx'];
%     nev_filename = [data_path 'RippleData\Q' date '\' monk(1) '_' date 'b_ICMS_' dimInfo '.nev'];
%     ns2_filename = [data_path 'RippleData\Q' date '\' monk(1) '_' date 'b_ICMS_' dimInfo '.ns2'];
%     h5_filename = [data_path 'UTCSData\Q' date '\' monk(1) '_' date 'b_ICMS_' dimInfo '.h5'];

    if exist(plx_filename, 'file')
        %Reading Plexon file (get markers and Xdsip and Ydisp Time series data)
        [marker_ts, marker_data] = FMLReadMarkers(plx_filename);
        ts_xDisp = FMLReadTimeseries(plx_filename, 'X_disp');
        ts_yDisp = FMLReadTimeseries(plx_filename, 'Y_disp');
    else
        %Reading Ripple file (get markers and Xdsip and Ydisp Time series data)
        try
            [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
        catch
            nev_filename = [data_path 'Ripple\' addPath monk(1) '_' date '_COT_ICMS_' dimInfo '_shrink.nev'];
            [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
        end
    %     marker_data = marker_data(1:3:end);
    %     marker_ts = marker_ts(1:3:end);
        try
            ts_xDisp = FMLReadTimeseries(ns2_filename, 'analog 7');
            ts_yDisp = FMLReadTimeseries(ns2_filename, 'analog 8');
        catch
            ts_xDisp = FMLReadTimeseries(ns2_filename, 'analog 3');
            ts_yDisp = FMLReadTimeseries(ns2_filename, 'analog 4');
        end
    end

    %Reading H5 file
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
    
end
    
if correctEventCodes
    %Remove duplicated event codes and zeros inserted after bit change
    nonBitChangeNDX = marker_data ~= 0; %Find non-zero markers
    nonEventCodes = [0; diff(nonBitChangeNDX)] == -1; %Find extra zeros after non-zero markers (extra zero only placed after non-zero bit change)
    dupEventCodes = marker_data(~nonEventCodes);
    temp1_marker_ts = marker_ts(~nonEventCodes);
    filterCodes = dupEventCodes(([1000; diff(dupEventCodes)] ~= 0 | dupEventCodes == 0)); %Remove non-zero duplicates only
    filterTS = temp1_marker_ts(([1000; diff(dupEventCodes)] ~= 0 | dupEventCodes == 0));

    %Reinsert codes for rare cases of true duplicates in true stream (251 251 may occur with trial start and serial ID # beginning with 251)
    trueCodes = cat(1,trial_data.event_codes);
    trueDup = find([1000; diff(trueCodes)] == 0);
    nnzDup = trueDup(trueCodes(trueDup)>0);
    for di = 1:length(nnzDup)
        filterCodes = [filterCodes(1:nnzDup(di)-1); trueCodes(nnzDup(di)); filterCodes(nnzDup(di):end)];
        tsNDX = find(marker_ts == filterTS(nnzDup(di)-1)) + 3; %Find data point of non-zero duplicate from marker_ts (3 points away for non-zer)
        filterTS = [filterTS(1:nnzDup(di)-1); marker_ts(tsNDX); filterTS(nnzDup(di):end)];
    end

    %Check to see stream was corrected properly
    if sum(filterCodes == trueCodes) ~= length(trueCodes)
        error('Event Code Stream not fixed properly')
    end

    marker_data = filterCodes;
    marker_ts = filterTS;

end
% trial_data_excluded = trial_data; % create trial_data_exclued variable -- protect trial_data variable from modifying

%% Target Extraction

%Get trial serial codes and instructed target from all trials
%Can get other information here as well if desired
numTrials = length(trial_data);
trialSerial = zeros(numTrials, 4);
trialTargets = zeros(numTrials, 1);
trialPriorSuccess = zeros(numTrials, 1);
trialKeep = zeros(numTrials, 1);
catchTrials = zeros(numTrials, 1);
blockNames = string(nan(numTrials, 1));

for i = 1:numTrials
    startCodeNDX = find(trial_data(i).event_codes == 251, 1, 'first');
    trialSerial(i, :) = trial_data(i).event_codes(startCodeNDX+1 : startCodeNDX+4);
    trialTargets(i) = trial_data(i).selected_target_index;
    trialPriorSuccess(i) = trial_data(i).event_codes(startCodeNDX+11);
    trialKeep(i) = sum(trialsUse == i);
    catchTrials(i) = trial_data(i).catch_trial;
    blockNames(i) =  string(trial_data(i).intstruction_block_name);
end

if isnan(trialsUse)
    trialKeep(:) = true;
end

%Get colors of each trial
[clrs] = COTReportTargetColorChange(h5_filename);

if matrixDim == "dim3"
    dimColor = "#87d387";
elseif matrixDim == "dim4"
    dimColor = "#9ecb9e";
elseif matrixDim == "dim5a"
    dimColor = "#b5c4b5";
elseif matrixDim == "dim5b"
    dimColor = "#b8c3b8";
elseif matrixDim == "dim5c"
    dimColor = "#bbc2bb";
elseif matrixDim == "dim5d"
    dimColor = "#bec1be";
elseif matrixDim == "dim5e"
    dimColor = "#bfc0bf";  
elseif matrixDim == "dimFull"
    dimColor = "#c0c0c0";
else
    dimColor = "All";
end

%Generate Markers for important patterns
startMarker = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.','.', '.', '.');
reachMarker = Marker(60);
holdMarker = Marker(248);
%errorMarker = Marker(100);
targetErrorMarker = Marker(100, 8);
holdErrorMarker = Marker(100, 9);

%Generate key patterns for extraction
holdCorrectPattern = Pattern(startMarker, reachMarker, holdMarker); %248 only provided if hold begins (correct target entered first)
holdErrorPattern = Pattern(startMarker, reachMarker, holdErrorMarker); %Correct target, left too soon
targetErrorPattern = Pattern(startMarker, reachMarker, targetErrorMarker); %Incorrect target

%Extract cursor position at marker times for each key pattern in X and Y
%Note that data is aligned on 'reach marker' and collected from 1 msec 
%before and after target indicator.  Data points were averaged 
%for determining target location
[T_correct_x, ~, PI_correct_x] = TimeseriesFindExtractAlign(ts_xDisp, marker_ts, marker_data, holdCorrectPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
[T_holdError_x, ~, PI_holdError_x] = TimeseriesFindExtractAlign(ts_xDisp, marker_ts, marker_data, holdErrorPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
[T_targetError_x, ~, PI_targetError_x] = TimeseriesFindExtractAlign(ts_xDisp, marker_ts, marker_data, targetErrorPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);

[T_correct_y, ~, PI_correct_y] = TimeseriesFindExtractAlign(ts_yDisp, marker_ts, marker_data, holdCorrectPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
[T_holdError_y, ~, PI_holdError_y] = TimeseriesFindExtractAlign(ts_yDisp, marker_ts, marker_data, holdErrorPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
[T_targetError_y, ~, PI_targetError_y] = TimeseriesFindExtractAlign(ts_yDisp, marker_ts, marker_data, targetErrorPattern,'AlignAtIndex',2,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);

%Assemble position and pattern data from each pattern type into single variable
alignData_x = {T_correct_x{1}, T_holdError_x{1}, T_targetError_x{1}};
alignData_y = {T_correct_y{1}, T_holdError_y{1}, T_targetError_y{1}};
pidx_x = {PI_correct_x{1}, PI_holdError_x{1}, PI_targetError_x{1}};
pidx_y = {PI_correct_y{1}, PI_holdError_y{1}, PI_targetError_y{1}};

%Variable initialization for loops
xData = [];
yData = [];
instructedTarget = [];

%Logic variables 
trialFound = [];     %For keeping track of trials not found (extra event codes at beginning)
isDimLevel = [];     %For removing trials with undesired dim levels
isBlock = [];           
isPriorSuccess = []; %For using only trials with prior success
isKeep = [];         %For specific ranges of trials
isCatch = [];        %For catch trials

% Combine x and y position of cursor at taget entry from every trial
for i = 1:length(alignData_x)
    for j = 1:length(alignData_x{i})
        
        %Get x and y data from current trial
        xData = [xData; mean(alignData_x{i}{j}.Data)];
        yData = [yData; mean(alignData_y{i}{j}.Data)];
        
        %Extract trial serial number from marker_data
        currentTrial = marker_data(pidx_x{i}(j,1)+1 : pidx_x{i}(j,1)+4)';
        
        %Correlate trial serial number with trial number from .h5 file (and clrs variable)
        trialNDX = sum(currentTrial == trialSerial, 2) == 4;
        
        %Keep track of trials non found
        trialFound = [trialFound; sum(trialNDX)];
%         if sum(trialNDX) == 0
%             badTrialCodes = [badTrialCodes; currentTrial'];
%         end
        
        %Call important variables from .h5 file (and clrs variable)
        instructedTarget = [instructedTarget; trialTargets(trialNDX)];
        isDimLevel = [isDimLevel; clrs(trialNDX) == dimColor];
        isBlock = [isBlock; strcmpi(blockNames(trialNDX), useBlock)];
        isPriorSuccess = [isPriorSuccess; trialPriorSuccess(trialNDX) == 31];
        isKeep = [isKeep; trialKeep(trialNDX)];
        isCatch = [isCatch; catchTrials(trialNDX)];
            
    end
end

if ~usePriorSuccess
    isPriorSuccess(:) = true;
    successStr = ' (Prior Errors Included)';
    saveStr_success = 'withPrevErrors';
else
    successStr = ' (Prior Success Only)';
    saveStr_success = 'prevSuccess';
end

%Reset dimColor if all levels are to be used
if dimColor == "All"
    isDimLevel(:) = 1;
end

%Skip trials not found in all necessary variables to make matrices same length
if sum(~trialFound > 0)
    warning([num2str(sum(~trialFound)) ' trial(s) could not be located by event codes and were ommitted'])
    xData(~trialFound) = [];
    yData(~trialFound) = [];
end

%Calculate angle of cursor at the moment a target was entered
enteredTarget = zeros(length(instructedTarget), 1);
cursorAngle = zeros(length(xData), 1);
for i = 1:length(xData)
    
    cursorAngle(i) = atan2(yData(i), xData(i));  
   
end

if numTargets == 8
    
    %Extrapolate angle into a target number
    enteredTarget(cursorAngle > -pi()/8 & cursorAngle <= pi()/8) = 0;
    enteredTarget(cursorAngle > pi()/8 & cursorAngle <= 3*pi()/8) = 1;
    enteredTarget(cursorAngle > 3*pi()/8 & cursorAngle <= 5*pi()/8) = 2;
    enteredTarget(cursorAngle > 5*pi()/8 & cursorAngle <= 7*pi()/8) = 3;
    enteredTarget(cursorAngle <= -7*pi()/8 | cursorAngle > 7*pi()/8) = 4;
    enteredTarget(cursorAngle <= -5*pi()/8 & cursorAngle > -7*pi()/8) = 5;
    enteredTarget(cursorAngle <= -3*pi()/8 & cursorAngle > -5*pi()/8) = 6;
    enteredTarget(cursorAngle <= -pi()/8 & cursorAngle > -3*pi()/8) = 7;

elseif numTargets == 4   
    enteredTarget(cursorAngle > -pi()/4 & cursorAngle <= pi()/4) = 0;
    enteredTarget(cursorAngle > pi()/4 & cursorAngle <= 3*pi()/4) = 1;
    enteredTarget(cursorAngle <= -3*pi()/4 | cursorAngle > 3*pi()/4) = 2;
    enteredTarget(cursorAngle <= -pi()/4 & cursorAngle > -3*pi()/4) = 3;
end


%Remove unwanted trials (can add other criteria besides dim level here in future)
% removeTrials = (~isDimLevel | ~isPriorSuccess | ~isKeep | ~isBlock);
removeTrials = (~isDimLevel | ~isPriorSuccess | ~isKeep);

catchEnter = enteredTarget(isCatch & ~removeTrials);
catchInstruct = instructedTarget(isCatch & ~removeTrials);

enteredTarget(removeTrials | isCatch) = [];
instructedTarget(removeTrials | isCatch) = [];
 

%% Matrix Generation
remLowTrialCount = false;
numLowTrials = 10;

%Count number of times each combination of instructed/entered targett occurred 
confusion = zeros(numTargets,numTargets);
catchConfuse = zeros(numTargets,numTargets);
for i = 1:numTargets
    for j = 1:numTargets
        
        confusion(i,j) = sum(instructedTarget == i-1 & enteredTarget == j-1);
        catchConfuse(i,j) = sum(catchInstruct == i-1 & catchEnter == j-1);
        
        
        sum(diag(confusion), 'all')/sum(confusion, 'all');
        sum(diag(catchConfuse), 'all')/sum(catchConfuse, 'all');

    end
end
trialsPerTarget = sum(confusion, 2);
trialsPerCatch = sum(catchConfuse, 2);

if remLowTrialCount
    if any(trialsPerTarget<=numLowTrials)
        lowNDX = trialsPerTarget<=numLowTrials
%         numTargets = numTargets - sum(lowNDX)
        confusion(lowNDX,:) = 0;
        trialsPerTarget (lowNDX) = 0;
    end
end

%Normalize confusion martix so sum(all entered targets) = 1 for each instructed target 
confuseStrength = zeros(size(confusion));
catchStrength = zeros(size(catchConfuse));
for i = 1:length(confusion)
    
    if sum(confusion(i,:)) == 0
        confuseStrength(i,:) = zeros(size(confusion(i,:)));  %Can't divide by zero
        catchStrength(i,:) = zeros(size(catchConfuse(i,:)));
    else
        confuseStrength(i,:) = confusion(i,:)/sum(confusion(i,:)); %Normalize to num of trials for each target
        catchStrength(i,:) = catchConfuse(i,:)/sum(catchConfuse(i,:));
    end
    
end



%% Plotting 
% close all

%Plot normalized confusion matrix with colored formatting
fontSize = 20;
showTitle = true;

figure('Position', [600 300 700 630])
for i = 1:numTargets
    for j = 1:numTargets   
        %patch([i-1 i i i-1], [j-1 j-1 j j], confuseStrength(i,j))
        cPtch(i,j) = patch([i-1 i i i-1], [j-1 j-1 j j], [1-confuseStrength(i,j), 1-confuseStrength(i,j), 1]);        
    end
end
step = 0.01;
colorChange = [0:step:1]';
colormap([1-colorChange 1-colorChange ones(size(colorChange))])
c = colorbar;
c.FontSize = fontSize;
xlim([0 numTargets])

% xticks([0.5:2:7.5])
% xticklabels({'T0','T2','T4','T6'}) %'T4','T5','T6','T7'})

if numTargets == 8   
    xticks([0.5 2.5 3.5 4.5 6.5 ])
    xticklabels({'T0','T2','T3','T4','T6'}) %'T4','T5','T6','T7'})
    yText = -0.75;
elseif numTargets == 4  
    xticks([0.5:numTargets])
    xticklabels({'T0','T1','T2','T3'})
    yText = -0.85;
end

xlabel({[''];['Instructed Target']}, 'FontSize', fontSize+16, 'fontWeight', 'Bold')
ylim([0 numTargets])
yticks([0.5:numTargets])
for yL = 1:numTargets
    yStrs{yL} = ['T' num2str(yL-1)];
end
yticklabels(yStrs) 
ylabel('Entered Target', 'FontSize', fontSize+18, 'fontWeight', 'Bold')
if showTitle
    tit = title(['Confusion Matrix, ' date ', ' titleStr successStr], 'Interpreter', 'none');
    tit.Position(2) = 8.2;
end
ax = gca;
ax.FontSize = fontSize;
ax.Title.FontSize = fontSize-6;
% ax.Position = [0.2 0.25 0.35 0.4];
% ax.YLabel.Position(1) = -0.9;
% ax.XLabel.Position(2) = -0.7;
% c.Position([1,2,4]) = [0.6 ax.Position([2,4])];
% c.Ticks = 0:0.2:1;

textNDX = sum(confusion,2) > 0;
% text(-0.9, yText, '(# Trials)', 'FontSize', fontSize);
textPos = [0.5:numTargets];
for i = 1:length(textNDX)
    if textNDX(i)
        text(textPos(i), yText, ['(' num2str(trialsPerTarget(i)) ')'], 'FontSize', fontSize-4, 'HorizontalAlignment', 'center')
    end
end

% BM(confusion)

saveas(gcf, [save_path 'confusionMat_' date '_' titleStr '_' saveStr_success '.png'])

round(255-(255*confuseStrength))'

%% Plot Catch
%Plot normalized catch trials with colored formatting
if sum(isCatch) > 0
    
    figure('Position', [600 300 700 630])
    for i = 1:numTargets
        for j = 1:numTargets
            catchPtch = patch([i-1 i i i-1], [j-1 j-1 j j], [1-catchStrength(i,j), 1-catchStrength(i,j), 1]);
        end
    end
    step = 0.01;
    colorChange = [0:step:1]';
    colormap([1-colorChange 1-colorChange ones(size(colorChange))])
    c = colorbar;
    c.FontSize = fontSize;
    xlim([0 numTargets])
    
    if numTargets == 8   
        xticks([0.5:2:numTargets ])
        xticklabels({'T0','T2','T4','T6'}) %'T4','T5','T6','T7'})
        yText = -0.75;
    elseif numTargets == 4  
        xticks([0.5:numTargets])
        xticklabels({'T0','T1','T2','T3'})    
        yText = -0.85;
    end

    xlabel({[''];['Instructed Target']}, 'FontSize', fontSize+16, 'fontWeight', 'Bold')
    ylim([0 numTargets])
    yticks([0.5:numTargets])
    for yL = 1:numTargets
        yStrs{yL} = ['T' num2str(yL-1)];
    end
    yticklabels(yStrs) 
    ylabel('Entered Target', 'FontSize', fontSize+18, 'fontWeight', 'Bold')
    if showTitle
        tit = title(['Catch Trials, ' date ', ' titleStr successStr], 'Interpreter', 'none');
        tit.Position(2) = 4.5;
    end
    ax = gca;
    ax.FontSize = fontSize;
    ax.Title.FontSize = fontSize-6;
    ax.Position = [0.2 0.25 0.35 0.4];
    ax.YLabel.Position(1) = -0.9;
    ax.XLabel.Position(2) = -0.7;
    c.Position([1,2,4]) = [0.6 ax.Position([2,4])];
    c.Ticks = [0:0.2:1];

    textNDX = sum(confusion,2) > 0;
    % text(-0.9, yText, '(# Trials)', 'FontSize', fontSize);
    textPos = [0.5:numTargets];
    for i = 1:length(textNDX)
        if textNDX(i)
            text(textPos(i), yText, ['(' num2str(trialsPerCatch(i)) ')'], 'FontSize', fontSize-4, 'HorizontalAlignment', 'center')
        end
    end

    saveas(gcf, [save_path 'confusionMat_' date '_' titleStr '_catchTrials_' saveStr_success '.png'])
    
end


