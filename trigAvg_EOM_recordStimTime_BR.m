close all; 
clear;

%% Initial Settings
addpath(genpath('D:\MATLAB\FML_matlab'))
addpath(genpath('D:\Projects\COT_ICMS\'))
addpath(genpath('D:\MATLAB\fileHandling_BR'))

monk = 'Qulio';
corticalRegion = 'AIP';

data_path = ['D:\RawData\COT_ICMS\' monk '\' corticalRegion '\'];
% data_path = ['D:\Projects\COT_ICMS\' monk '\' corticalRegion '\stimTA\emgData\'];
save_path = ['D:\Projects\COT_ICMS\' monk '\' corticalRegion '\stimTA\'];

EOM_date = '20211210';
ECal_date = '20211210';
fileToken = [monk(1) '_' EOM_date '_COT_ICMS_M1catch_78'];
artPath = ['D:\Projects\COT_ICMS\' monk '\' corticalRegion '\stimTA\emgData\'];
artToken = ['Q_20211210_COT_ICMS_M1catch_78'];

useNs2 = true; %1000 Hz
usePlx = false;

recArtTs = true;
artChan = 81;
useMC = true; %Use marker code for instruction onset rather than artifact times

trialNorm = 0;
trialRange = [1 10000] - trialNorm;
stimAnnulus = NaN;

% task_targets = [0, 2, 4, 6];
task_targets = [0, 1, 2, 3];

snipRange_EOM = [-200 800]/1000; % In sec


%% Data Extraction
%Reading H5 file

trialRange(trialRange<1) = 1;
useTrials = [trialRange(1):trialRange(2)];

if ~usePlx
    h5_filename = [data_path 'UTCS\' fileToken '.h5'];
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
    
    if max(useTrials) > length(trial_data)
        useTrials = useTrials(1):length(trial_data);
    end
    trialAnnuli = [trial_data.selected_annulus_index]';
    
    keepTrials = false(length(trial_data),1);
    keepTrials(useTrials) = true;
    if ~isnan(stimAnnulus)
        keepAnnuli = trialAnnuli == stimAnnulus;
    else
        keepAnnuli = true(size(trialAnnuli));
    end
    
    allTrials = keepAnnuli & keepTrials;
    h5_TC = {trial_data(allTrials).event_codes}';
    useTrial_SN = cellfun(@(x) x(find(x==251, 1, 'first')+1:find(x==251, 1, 'first')+4)', h5_TC, 'UniformOutput', false);
    useTrial_SN = cat(1, useTrial_SN{:});
end

if usePlx
    %Reading Plexon Files
    plx_filename = [data_path fileToken '.plx'];
    [marker_ts, marker_data] = FMLReadMarkers(plx_filename);
    EOM_fullFile = [data_path fileToken '.plx'];
else
    %Reading Ripple files
    try
        nev_filename = [data_path 'Ripple\' fileToken '.nev'];
        [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
        EOM_fullFile = [data_path 'Ripple\' fileToken '.ns2'];
    catch
        try
            nev_filename = [data_path 'Ripple\' fileToken '_shrink.nev'];
            [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
            EOM_fullFile = [data_path 'Ripple\' fileToken '.ns2'];
        catch
            try
                nev_filename = [data_path 'Ripple\Extra\' fileToken '.nev'];
                [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
                EOM_fullFile = [data_path 'Ripple\Extra\' fileToken '.ns2'];
            catch
                nev_filename = [data_path 'Ripple\Extra\' fileToken '_shrink.nev'];
                [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
                EOM_fullFile = [data_path 'Ripple\Extra\' fileToken '.ns2'];
            end
        end
    end
end

if ~usePlx

    if recArtTs
        sma_chan = 2;
        artTs = FMLReadNEVDIOEventsTimestamps(nev_filename, sma_chan);
        if isempty(artTs)
            artTs = FMLReadNEVDIOEventsTimestamps(nev_filename, sma_chan+1);
        end
        if isempty(artTs)
            error('Need artifact file - no stim times recorded')
        end
        ISI = diff(artTs);
        repeatNDX = [false ISI<0.001];
        artTs(repeatNDX) = [];
    elseif ~useMC
        art_filename = [artPath artToken '.nex'];
        nexArtInfo = readNexFile(art_filename);
        nexArtInfo.channels = cellfun(@(x) sscanf(x.name,'sig%d'), nexArtInfo.neurons);
        nexArtInfo.units = cellfun(@(x) x.unitNumber, nexArtInfo.neurons);
        nexArtInfo.tscounts = cellfun(@(x) length(x.timestamps), nexArtInfo.neurons);
        nexArtInfo.wfcounts = cellfun(@(x) size(x.waveforms,2), nexArtInfo.waves);
    
        keepNDX = cellfun(@(x) ~strcmpi(x.name(end), 'U'), nexArtInfo.neurons);
    
        nexArtInfo.channels(~keepNDX) = [];
        nexArtInfo.units(~keepNDX) = [];
        nexArtInfo.neurons(~keepNDX) = [];
        nexArtInfo.waves(~keepNDX) = [];
        nexArtInfo.wfcounts(~keepNDX) = [];
        nexArtInfo.tscounts(~keepNDX) = [];
    
        artStruct = cat(1,nexArtInfo.neurons{:});
        artTs = sort(cat(1,artStruct(:).timestamps));
    end
    
    tsEOM_x = FMLReadTimeseries(EOM_fullFile, 'analog 5');
    EOM_data(1,:) = tsEOM_x.Data;
    tsEOM_y = FMLReadTimeseries(EOM_fullFile, 'analog 6');
    EOM_data(2,:) = tsEOM_y.Data;
    ts = tsEOM_x.StartTime;
    adfreq = 1/tsEOM_x.Period;

else    
    %Reading Plexon file (get markers and Xdsip and Ydisp Time series data)
    art_filename = [artPath artToken '.plx'];
    [~, ~, artTs_1, ~] = plx_waves_v(art_filename, artChan, 1);
    [~, ~, artTs_2, ~] = plx_waves_v(art_filename, artChan, 2);
    [~, ~, artTs_3, ~] = plx_waves_v(art_filename, artChan, 3);
    [~, ~, artTs_4, ~] = plx_waves_v(art_filename, artChan, 4);   
    artTs = [artTs_1; artTs_2; artTs_3; artTs_4];

    [adfreq, ~, ts, ~, EOM_data(1,:)] = plx_ad_v(EOM_fullFile, 'EOM_H');
    [~, ~, ~, ~, EOM_data(2,:)] = plx_ad_v(EOM_fullFile, 'EOM_V');

    tsEOM_x = UniformTimeseries();
    tsEOM_x.Period = 1/adfreq;
    tsEOM_x.Data = EOM_data(1,:);
    tsEOM_x.StartTime = ts;

    tsEOM_y = UniformTimeseries();
    tsEOM_y.Period = 1/adfreq;
    tsEOM_y.Data = EOM_data(2,:);
    tsEOM_y.StartTime = ts;

end



%% Triggered Averaging

tAnalog = [1/adfreq:1/adfreq:length(EOM_data)/adfreq]' + ts;

allTargets = task_targets;
startMarker = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.');
instructMarker = Marker(210);
holdMarker = Marker(248);
endMarker = Marker(253);

for targ_i = 1:length(allTargets)

    targetMarker = Marker(allTargets(targ_i));
    targetPattern = Pattern(startMarker, targetMarker, instructMarker, holdMarker, endMarker);
    patNDX = FindPattern(marker_data, targetPattern);

    all_EOM_snips_stimOn = [];
    all_EOM_snips_stimOff = [];
    if useMC
        
        [T_EOM_xOn, ~, PI_EOM_xOn] = TimeseriesFindExtractAlign(tsEOM_x, marker_ts, marker_data, targetPattern,'AlignAtIndex',3,'ExtractMode','around',...
                                     'ExtractAroundIndex',3,'ExtractBeforeDuration',abs(snipRange_EOM(1)), 'ExtractAfterDuration', abs(snipRange_EOM(2))+0.001);
        [T_EOM_yOn, ~, PI_EOM_yOn] = TimeseriesFindExtractAlign(tsEOM_y, marker_ts, marker_data, targetPattern,'AlignAtIndex',3,'ExtractMode','around',...
                                     'ExtractAroundIndex',3,'ExtractBeforeDuration',abs(snipRange_EOM(1)), 'ExtractAfterDuration', abs(snipRange_EOM(2))+0.001);
        
        snipSize = length(snipRange_EOM(1):1/adfreq:snipRange_EOM(2));
        xySnips_On = cellfun(@(x,y) [x.Data(1:snipSize); y.Data(1:snipSize)], T_EOM_xOn{1}(:), T_EOM_yOn{1}(:), 'UniformOutput', false);        
        all_EOM_snips_stimOn = cat(3, xySnips_On{:});
        
        [T_EOM_xOff, ~, PI_EOM_xOff] = TimeseriesFindExtractAlign(tsEOM_x, marker_ts, marker_data, targetPattern,'AlignAtIndex',4,'ExtractMode','around',...
                                     'ExtractAroundIndex',4,'ExtractBeforeDuration',abs(snipRange_EOM(1)), 'ExtractAfterDuration', abs(snipRange_EOM(2))+0.001);
        [T_EOM_yOff, ~, PI_EOM_yOff] = TimeseriesFindExtractAlign(tsEOM_y, marker_ts, marker_data, targetPattern,'AlignAtIndex',4,'ExtractMode','around',...
                                     'ExtractAroundIndex',4,'ExtractBeforeDuration',abs(snipRange_EOM(1)), 'ExtractAfterDuration', abs(snipRange_EOM(2))+0.001);
        
        xySnips_Off = cellfun(@(x,y) [x.Data(1:snipSize); y.Data(1:snipSize)], T_EOM_xOff{1}(:), T_EOM_yOff{1}(:), 'UniformOutput', false);        
        all_EOM_snips_stimOff = cat(3, xySnips_Off{:});
    
    else 
        
        stimDur = cell(length(allTargets),1);
        for ti = 1:size(patNDX, 1)

           trialTime = [marker_ts(patNDX(ti,1)) marker_ts(patNDX(ti,end))];
           art_ts_ndx = (artTs >= trialTime(1) & artTs <= trialTime(2));
           trialArt_ts = artTs(art_ts_ndx);

           snipStart1_EOM = find(tAnalog >= min(trialArt_ts) + snipRange_EOM(1), 1, 'first');
           snipEnd1_EOM = snipStart1_EOM + (snipRange_EOM(2)-snipRange_EOM(1)) * adfreq ;
           EOM_snip1 = EOM_data(:, snipStart1_EOM : snipEnd1_EOM);
           all_EOM_snips_stimOn = cat(3, all_EOM_snips_stimOn, EOM_snip1);
           
           snipStart2_EOM = find(tAnalog >= max(trialArt_ts) + snipRange_EOM(1), 1, 'first');
           snipEnd2_EOM = snipStart2_EOM + (snipRange_EOM(2)-snipRange_EOM(1)) * adfreq ;
           EOM_snip2 = EOM_data(:, snipStart2_EOM : snipEnd2_EOM);
           all_EOM_snips_stimOff = cat(3, all_EOM_snips_stimOff, EOM_snip2);
           
           stimDur{targ_i}(ti,1) = trialArt_ts(end)-trialArt_ts(1);
           
        end 
        
    end
    
    if ~usePlx
        nev_SN = [marker_data(patNDX(:,1)+1), marker_data(patNDX(:,1)+2), marker_data(patNDX(:,1)+3), marker_data(patNDX(:,1)+4)];
        snMatch = false(size(patNDX,1),1);
        for pi = 1:size(patNDX,1)
            snMatch(pi,1) = any(sum(nev_SN(pi,:) == useTrial_SN,2)==4);
        end
        all_EOM_snips_stimOn(:,:,~snMatch) = [];
        all_EOM_snips_stimOff(:,:,~snMatch) = [];
        
        if exist('stimDur', 'var')
            stimDur{targ_i}(~snMatch, :) = [];
        end
    end
    
    snipsEOM_On{targ_i} = all_EOM_snips_stimOn;
    snipsEOM_Off{targ_i} = all_EOM_snips_stimOff;
    nTrials(targ_i) = size(all_EOM_snips_stimOn, 3);
    

%     figure
%     plot(mean(all_artWaves, 1));

%     yline(2*std(stimTA) + mean(stimTA), 'r');
%     yline(-(2*std(stimTA)) + mean(stimTA), 'r');
    
end

%% Eye Movement Calibration

EOM_data_path = ['D:\RawData\COT_ICMS\' monk '\' corticalRegion '\'];
EOM_save_path = ['D:\Projects\COT_ICMS\' monk '\' corticalRegion '\StimTA\'];

EOM_ns2_file = [EOM_data_path 'Ripple\' monk(1) '_' ECal_date '_eyemvmt.ns2'];
EOM_nev_file = [EOM_data_path 'Ripple\' monk(1) '_' ECal_date '_eyemvmt.nev'];
EOM_plx_file = [EOM_data_path 'Plexon\' monk(1) '_' ECal_date '_eyemvmt.plx'];

if ~usePlx
    [marker_ts2, marker_data2] = FMLReadMarkers(EOM_nev_file);
    tsEOM_x2 = FMLReadTimeseries(EOM_ns2_file, 'analog 5');
    tsEOM_y2 = FMLReadTimeseries(EOM_ns2_file, 'analog 6');
else
    [marker_ts2, marker_data2] = FMLReadMarkers(EOM_plx_file);
    [adfreq, ~, ts, ~, tsEOM_x2] = plx_ad_v(EOM_plx_file, 'EOM_H');
    [~, ~, ~, ~, tsEOM_y2]= plx_ad_v(EOM_plx_file, 'EOM_V');
end

eyeTargets = [0, 3, 6, 9];
% eyeTargets = [0, 4, 8, 12];
eyePos = zeros(length(eyeTargets),2);

for ti = 1:length(eyeTargets)
    
    startTrial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.');
%     successMarker = Marker(119);
    enterMarker = Marker(248);
    targetMarker2 = Marker(eyeTargets(ti));
    targetPattern2 = Pattern(startTrial, targetMarker2, enterMarker);%, successMarker);
    
    xCollect = [];
    yCollect = [];

    if ~usePlx

        [T_x, ~, PI_x] = TimeseriesFindExtractAlign(tsEOM_x2, marker_ts2, marker_data2, targetPattern2,'AlignAtIndex',1,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
        [T_y, ~, PI_y] = TimeseriesFindExtractAlign(tsEOM_y2, marker_ts2, marker_data2, targetPattern2,'AlignAtIndex',1,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
    
        for utsi = 1:length(T_x{1})
            
            xCollect = [xCollect; T_x{1}{utsi}.Data'];
            yCollect = [yCollect; T_y{1}{utsi}.Data'];
            
        end

    else
        
        time = [1/adfreq : 1/adfreq : length(tsEOM_x2)/adfreq]+ts;
        patNDX2 = FindPattern(marker_data2, targetPattern2);
        timeNDX = zeros(size(patNDX2,1),1);
        for time_i = 1:length(timeNDX)
            normTime = abs(time-marker_ts2(patNDX2(time_i,3)));
            timeNDX = find(normTime == min(normTime))+1;
            
            xCollect = [xCollect; tsEOM_x2(timeNDX)];
            yCollect = [yCollect; tsEOM_y2(timeNDX)];
        end

    end

        eyePos(ti, 1) = mean(xCollect);
        eyePos(ti, 2) = mean(yCollect); 

end

% xOffset = mean(eyePos([2,4],1));
% yOffset = mean(eyePos([1,3],2));
% eyePos = [eyePos(:,1)-xOffset, eyePos(:,2)-yOffset];
eyePos(isnan(eyePos(:,1)),:) = [];

% annuliDist = 2.875;
annuliDist = 2.125;
screenDist = 36.25;
% screenDist = 48;
numTargets = 16;
targetNums = [0:4:15];

targetVolt = mean(max(abs(eyePos),[],2));
degPerVolt = atand(annuliDist/screenDist) / targetVolt;

scatter(eyePos(:,1)*degPerVolt, eyePos(:,2)*degPerVolt, 'k', 'filled', 'SizeData', 100);
hold on
for ti = 1:length(targetNums)
    [xPtch, yPtch] = getTargetPatch(annuliDist, annuliDist+1.5, numTargets, targetNums(ti));
    patch(xPtch, yPtch, 'k', 'FaceAlpha', 0.2)
end
hold off
axis([-5 5 -5 5])

if exist('stimDur', 'var')
    save([save_path 'EOM_snipInfo_' EOM_date '_' num2str(adfreq) 'Hz.mat'], 'snipsEOM_On', 'snipsEOM_Off', 'stimDur', 'allTargets', 'snipRange_EOM', 'adfreq', 'degPerVolt', '-v7.3')
else
    save([save_path 'EOM_snipInfo_' EOM_date '_' num2str(adfreq) 'Hz.mat'], 'snipsEOM_On', 'snipsEOM_Off', 'allTargets', 'snipRange_EOM', 'adfreq', 'degPerVolt', '-v7.3')
end

 %% Plotting
close all
adfreq = 1000;
load([save_path 'EOM_snipInfo_' EOM_date '_' num2str(adfreq) 'Hz.mat']);

% monk = 'Felix';
% corticalRegion = 'S1';
% matFile = 'EOM_snipInfo_20211005_1000Hz.mat';
% matDir = ['D:\Projects\COT_ICMS\' monk '\' corticalRegion '\stimTA\'];
% load([matDir matFile])

pltRange_EOM = [-200 500];

EOM_onFig = figure('Position', [500 300 900 600]);
EOM_offFig = figure('Position', [500 300 900 600]);

axFS= 10;
axLW = 2;

eyeNames = ["X Poistion (deg)", "Y Poistion (deg)"];

plotTime_EOM = [snipRange_EOM(1) : 1/adfreq : snipRange_EOM(2)]*1000;
plotTime_NDX = plotTime_EOM >= pltRange_EOM(1) & plotTime_EOM <= pltRange_EOM(2);
plotTime_EOM = plotTime_EOM(plotTime_NDX);

ylim_EOM = [-20 20];
xlim_EOM = [pltRange_EOM(1)-5, pltRange_EOM(2)+5];

switch monk(1)
    case 'Q'
        centerRad = 3.54/2;
        innerRad = 2.76;
        outerRad = innerRad+2.76;
        numTargets = 8;
        targetNums = [0 2 4 6];
    case 'F'
        centerRad = 3.27/2;
        innerRad = 2.76;
        outerRad = innerRad+2.36;
        numTargets = 4;
        targetNums = [0 1 2 3];
end


for targ_i = 1:length(allTargets)

    if ~isempty(snipsEOM_On{targ_i})
    
        stimOn_Trig_EOM = (snipsEOM_On{targ_i});
        stimOn_Trig_EOM = stimOn_Trig_EOM(:, plotTime_NDX, :);
        stimOff_Trig_EOM = (snipsEOM_Off{targ_i});
        stimOff_Trig_EOM = stimOff_Trig_EOM(:, plotTime_NDX, :);
    
        figure(EOM_onFig)
        onPlot_x = reshape(stimOn_Trig_EOM(1,:,:), size(stimOn_Trig_EOM(1,:,:), 2,3))*degPerVolt;
        onPlot_y = reshape(stimOn_Trig_EOM(2,:,:), size(stimOn_Trig_EOM(2,:,:), 2,3))*degPerVolt;
        xSat = sum(onPlot_x < (min(onPlot_x, [], 'all')+1) ,1) > 0;
        ySat = sum(onPlot_y < (min(onPlot_y, [], 'all')+1) ,1) > 0;
        satNDX = xSat & ySat;
    
%         color = zeros(length(plotTime_EOM), 3);
%         clrTrans = 0: 1/(round(length(color)/2)-1) :1; 
%         
%         color(1:length(clrTrans), 1) = flip(clrTrans);
%         color(1:length(clrTrans), 2) = 1;
%         
%     %     color(length(clrTrans):end,1) = (clrTrans);
%         color(length(clrTrans):end,2) = flip(clrTrans)/4 + 0.75;
%         color(length(clrTrans):end,2) = flip(clrTrans);
%         color(length(clrTrans):end,3) = clrTrans;
%     
%     %     color = flip(color); 
%         color(:,2) = color(:,2); color(color>1) = 1; color(color<0) = 0;
%         color = flip(color);
%         
%     %     color(plotTime_EOM==0, :) = [0 0 0];
%     %     color(plotTime_EOM > 0, 1) = 1;
        cMap = parula;
%         color = linspace(1,10,length(plotTime_EOM));
        color = 'k';

        EOM_Axes_on{1, targ_i} = subplot(3,length(allTargets),targ_i);
        hold on
        EOM_Plot_on{1, targ_i} = plot(plotTime_EOM, onPlot_x(:,~satNDX), 'color', [0,0,0]+0.8);
%         xline(0, '--r', 'ICMS On', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
        xline(0, '--r', 'LineWidth', 1.5);
        meanEOM_x = mean(onPlot_x, 2);
        scatter(plotTime_EOM, meanEOM_x, 5, color, 'filled');
%         colormap(cMap)
        ylim(ylim_EOM)
        xlim(xlim_EOM)
        if targ_i == 1
            ylabel("X Poistion (deg)")
        end
        title(['Target ' num2str(allTargets(targ_i)) ' (n= ' num2str(size(snipsEOM_On{targ_i}, 3)-sum(satNDX)) ')'])
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
    
        EOM_Axes_on{2, targ_i} = subplot(3,length(allTargets),targ_i+length(allTargets));
        hold on
        EOM_Plot_on{2, targ_i} = plot(plotTime_EOM, onPlot_y(:,~satNDX), 'color', [0,0,0]+0.8);
%         xline(0, '--r', 'ICMS On', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
        xline(0, '--r', 'LineWidth', 1.5);
        meanEOM_y = mean(onPlot_y, 2);
        scatter(plotTime_EOM, meanEOM_y, 5, color, 'filled');
%         colormap(cMap)
        hold off
        ylim(ylim_EOM)
        xlim(xlim_EOM)
        xlabel('Time (msec)')
        if targ_i == 1
            ylabel("Y Poistion (deg)")
        end
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
    
        EOM_Axes_on{3,targ_i} = subplot(3,length(allTargets),2*length(allTargets)+targ_i);
        hold on
        EOM_Plot_on{3,targ_i} = plot(onPlot_x(:,~satNDX), onPlot_y(:,~satNDX), 'color', [0,0,0]+0.9);
    %     scatter(meanEOM_x, meanEOM_y, 0.001*linspace(1,10*length(meanEOM_x),length(meanEOM_x)), color, 'filled');
        scatter(meanEOM_x, meanEOM_y, 5, color, 'filled');
%         colormap(cMap)
        [xPtch, yPtch] = getTargetPatch(innerRad, outerRad, numTargets, targetNums(targ_i));
        patch(xPtch, yPtch, 'k', 'FaceAlpha', 0.1);
        [xPtch2, yPtch2] = getTargetPatch(0, centerRad, 1, 0);
        patch(xPtch2(2:end), yPtch2(2:end), 'k', 'FaceAlpha', 0.1);
        hold off
        xlim(ylim_EOM)
        ylim(ylim_EOM)
        xlabel('X Position (deg)')
        if targ_i == 1
            ylabel('Y Position (deg)')
        end
        axis square
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
       
        figure(EOM_offFig)
        offPlot_x = reshape(stimOff_Trig_EOM(1,:,:), size(stimOff_Trig_EOM(1,:,:), 2,3))*degPerVolt;
        offPlot_y = reshape(stimOff_Trig_EOM(2,:,:), size(stimOff_Trig_EOM(2,:,:), 2,3))*degPerVolt;
        xSat = sum(offPlot_x < (min(offPlot_x, [], 'all')+1) ,1) > 0;
        ySat = sum(offPlot_y < (min(offPlot_y, [], 'all')+1) ,1) > 0;
        satNDX = xSat & ySat;
        
        EOM_Axes_off{1, targ_i} = subplot(3,length(allTargets),targ_i);
        hold on
        EOM_Plot_off{1, targ_i} = plot(plotTime_EOM, offPlot_x(:,~satNDX), 'color', [0,0,0]+0.9);
        xline(0, '--r', 'ICMS Off', 'LineWidth', 1, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
        meanEOM_x = mean(offPlot_x, 2);
        scatter(plotTime_EOM, meanEOM_x, 5, color, 'filled');
        hold off
        ylim(ylim_EOM)
        xlim(xlim_EOM)
        if targ_i == 1
            ylabel("X Poistion (deg)")
        end
        title(['Target ' num2str(allTargets(targ_i)) ' (n= ' num2str(size(snipsEOM_On{targ_i}, 3)-sum(satNDX)) ')'])
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
    
        EOM_Axes_off{2, targ_i} = subplot(3,length(allTargets),targ_i+length(allTargets));
        hold on
        EOM_Plot_off{2, targ_i} = plot(plotTime_EOM, offPlot_y(:,~satNDX), 'color', [0,0,0]+0.9);
        xline(0, '--r', 'ICMS Off', 'LineWidth', 1, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
        meanEOM_y = mean(offPlot_y, 2);
        scatter(plotTime_EOM, meanEOM_y, 5, color, 'filled');
        hold off
        ylim(ylim_EOM)
        xlim(xlim_EOM)
        xlabel('Time (msec)')
        if targ_i == 1
            ylabel("Y Poistion (deg)")
        end
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
    
        EOM_Axes_off{3,targ_i} = subplot(3,length(allTargets),2*length(allTargets)+targ_i);
        hold on
        EOM_Plot_off{3,targ_i} = plot(offPlot_x(:,~satNDX), offPlot_y(:,~satNDX), 'color', [0,0,0]+0.9);
    %     scatter(meanEOM_x, meanEOM_y, 0.001*linspace(1,10*length(meanEOM_x),length(meanEOM_x)), color, 'filled');
        scatter(meanEOM_x, meanEOM_y, 5, color, 'filled');
        [xPtch, yPtch] = getTargetPatch(innerRad, outerRad, numTargets, targetNums(targ_i));
        patch(xPtch, yPtch, 'k', 'FaceAlpha', 0.1)
        [xPtch2, yPtch2] = getTargetPatch(0, centerRad, 1, 0);
        patch(xPtch2(2:end), yPtch2(2:end), 'k', 'FaceAlpha', 0.1);
        hold off
        xlim(ylim_EOM)
        ylim(ylim_EOM)
        xlabel('X Position (deg)')
        if targ_i == 1
            ylabel('Y Position (deg)')
        end
        axis square
        box on
        ax = gca;
        ax.FontSize = axFS;
        ax.LineWidth = axLW;
    end  
end

figure(EOM_onFig)
% sgtitle('Stimulus-Triggered Averages of Eye Movements at ICMS Onset');
EOM_Axes_on{1,1}.Position([1,2,3]) = [0.05, 0.74, 0.2]; EOM_Axes_on{2,1}.Position([1,2,3]) = [0.05, 0.45, 0.2]; EOM_Axes_on{3,1}.Position = [0, 0.07, 0.3, 0.3];
EOM_Axes_on{1,2}.Position([1,2,3]) = [0.29, 0.74, 0.2]; EOM_Axes_on{2,2}.Position([1,2,3]) = [0.29, 0.45, 0.2]; EOM_Axes_on{3,2}.Position = [0.24, 0.07, 0.3, 0.3];
EOM_Axes_on{1,3}.Position([1,2,3]) = [0.53, 0.74, 0.2]; EOM_Axes_on{2,3}.Position([1,2,3]) = [0.53, 0.45, 0.2]; EOM_Axes_on{3,3}.Position = [0.48, 0.07, 0.3, 0.3];
EOM_Axes_on{1,4}.Position([1,2,3]) = [0.77, 0.74, 0.2]; EOM_Axes_on{2,4}.Position([1,2,3]) = [0.77, 0.45, 0.2]; EOM_Axes_on{3,4}.Position = [0.72, 0.07, 0.3, 0.3];
print([save_path monk(1) '_' EOM_date '_EOMstimTA_Onset_annulus' num2str(stimAnnulus) '.tiff'], '-dtiff', '-r300')

figure(EOM_offFig)
% sgtitle('Stimulus-Triggered Averages of Eye Movements at ICMS Offset');
EOM_Axes_off{1,1}.Position([1,2,3]) = [0.05, 0.74, 0.2]; EOM_Axes_off{2,1}.Position([1,2,3]) = [0.05, 0.45, 0.2]; EOM_Axes_off{3,1}.Position = [0, 0.07, 0.3, 0.3];
EOM_Axes_off{1,2}.Position([1,2,3]) = [0.29, 0.74, 0.2]; EOM_Axes_off{2,2}.Position([1,2,3]) = [0.29, 0.45, 0.2]; EOM_Axes_off{3,2}.Position = [0.24, 0.07, 0.3, 0.3];
EOM_Axes_off{1,3}.Position([1,2,3]) = [0.53, 0.74, 0.2]; EOM_Axes_off{2,3}.Position([1,2,3]) = [0.53, 0.45, 0.2]; EOM_Axes_off{3,3}.Position = [0.48, 0.07, 0.3, 0.3];
EOM_Axes_off{1,4}.Position([1,2,3]) = [0.77, 0.74, 0.2]; EOM_Axes_off{2,4}.Position([1,2,3]) = [0.77, 0.45, 0.2]; EOM_Axes_off{3,4}.Position = [0.72, 0.07, 0.3, 0.3];

saveas(EOM_onFig, [save_path monk(1) '_' EOM_date '_EOMstimTA_Onset_annulus' num2str(stimAnnulus) '.png'])

% saveas(EOM_onFig, [save_path monk(1) '_' EOM_date '_EOMstimTA_Onset_trials' num2str(trialRange(1)), '_' num2str(trialRange(2)) '.fig']) 
saveas(EOM_offFig, [save_path monk(1) '_' EOM_date '_EOMstimTA_Offset_annulus' num2str(stimAnnulus) '.png'])
% saveas(EOM_offFig, [save_path monk(1) '_' EOM_date '_EOMstimTA_Offset_trials' num2str(trialRange(1)), '_' num2str(trialRange(2)) '.fig']) 

 
function [xPtch, yPtch] = getTargetPatch(innerRad, outerRad, numTargets, targetNum)

targetSize = 360/numTargets;
theta = [-targetSize/2 : 0.01 : targetSize/2] + (targetSize*targetNum);

if innerRad == 0 && numTargets == 1
    xPtch = outerRad*cosd(theta);
    yPtch = outerRad*sind(theta);   
else
    xInner = innerRad*cosd(theta);
    yInner = innerRad*sind(theta);
    
    xOuter = outerRad*cosd(theta);
    yOuter = outerRad*sind(theta);
    
    xPtch = [xInner, flip(xOuter)]';
    yPtch = [yInner, flip(yOuter)]';
end


end
