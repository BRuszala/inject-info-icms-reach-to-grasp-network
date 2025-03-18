close all; clear;

%% Initial Settings
addpath(genpath('F:\MATLAB\FML_matlab'))
addpath(genpath('F:\MATLAB\fileHandling_BR'))

monk = 'Qulio';
cortical_region = 'PMv';

data_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\stimTA\emgData\'];
save_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\StimTA\'];

date = '20210414';

useNs2 = true;
useNs5 = false;
ns2_filename = 'Q_20211210_COT_ICMS_M1catch_78.ns2';
ns5_filename = 'Q_20211210_COT_ICMS_M1catch_78.ns5';
plx_filename = 'Q_20210415_ICMS_Dim4_EMG.plx';

task_targets = [0, 2, 4, 6];
EMG_target = [];
M1_elec = '001';

snipRange_EMG = [-30 50]/1000; % In sec
snipRange_EOM = [-200 800]/1000; % In sec

art_data_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\stimTA\emgData\'];
artChan = [001]; % For plexon files
art_filename = 'Q_20210414_ICMS_Dim4_EMG_art001.plx';
% art_filename = 'Q_20210416_ICMS_Dim4_EMG_art01.plx';
% filename = 'Q_20210326_ABCTEFYZ_EMG_Autosort.nex';

%% Data Extraction

art_fullFile = [art_data_path art_filename];
[marker_ts, marker_data] = FMLReadMarkers(art_fullFile);

%Reading Plexon file (get markers and Xdsip and Ydisp Time series data)
[~, ~, artTs_1, wave_1] = plx_waves_v(art_fullFile, ['sig' num2str(artChan, '%03.0f')], 1);
[~, ~, artTs_2, wave_2] = plx_waves_v(art_fullFile, ['sig' num2str(artChan, '%03.0f')], 2);
[~, ~, artTs_3, wave_3] = plx_waves_v(art_fullFile, ['sig' num2str(artChan, '%03.0f')], 3);
[~, ~, artTs_4, wave_4] = plx_waves_v(art_fullFile, ['sig' num2str(artChan, '%03.0f')], 4);
    
artTs = [artTs_1; artTs_2; artTs_3; artTs_4];
artWaves = [wave_1; wave_2; wave_3; wave_4];

plx_fullFile = [data_path plx_filename];
[adfreq, n, ts, fn, EMG_data(1,:)] = plx_ad_v(plx_fullFile, 'AD09');
[~, ~, ~, ~, EMG_data(2,:)] = plx_ad_v(plx_fullFile, 'AD10');
[~, ~, ~, ~, EMG_data(3,:)] = plx_ad_v(plx_fullFile, 'AD11');
[~, ~, ~, ~, EMG_data(4,:)] = plx_ad_v(plx_fullFile, 'AD12');
[~, ~, ~, ~, EMG_data(5,:)] = plx_ad_v(plx_fullFile, 'AD13');
[~, ~, ~, ~, EMG_data(6,:)] = plx_ad_v(plx_fullFile, 'AD14');

[~, ~, ~, ~, EOM_data(1,:)] = plx_ad_v(plx_fullFile, 'X_raw');
[~, ~, ~, ~, EOM_data(2,:)] = plx_ad_v(plx_fullFile, 'Y_raw');

% rectAD = abs(ad);
% smooth_rectAD = conv();
EMG_data = abs(EMG_data);
%Add smoothing?

%% Triggered Averaging

tAnalog = [1/adfreq:1/adfreq:length(EMG_data)/adfreq]' + ts;

allTargets = [task_targets EMG_target];
startMarker = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.');
instructMarker = Marker(210);
endMarker = Marker(253);

for ii = 1:length(allTargets)

    targetMarker = Marker(allTargets(ii));
    targetPattern = Pattern(startMarker, targetMarker, instructMarker, endMarker);

    patNDX = FindPattern(marker_data, targetPattern);

    all_EMG_snips = [];
    all_EOM_snips = [];
    all_artWaves = [];
    for ti = 1:size(patNDX, 1)

       trialTime = [marker_ts(patNDX(ti,1)) marker_ts(patNDX(ti,end))];
       art_ts_ndx = (artTs >= trialTime(1) & artTs <= trialTime(2));
       trialArt_ts = artTs(art_ts_ndx);
       all_artWaves = [all_artWaves; artWaves(art_ts_ndx, :)];
       
       for ai = 1:length(trialArt_ts)

           snipStart_EMG = find(tAnalog >= trialArt_ts(ai)+snipRange_EMG(1), 1, 'first');
           snipEnd_EMG = snipStart_EMG + (snipRange_EMG(2)-snipRange_EMG(1)) * adfreq;
           
           EMG_snip = EMG_data(:, snipStart_EMG : snipEnd_EMG);
           all_EMG_snips = cat(3, all_EMG_snips, EMG_snip);
           
       end
       
       if ~isempty(trialArt_ts)
           
           snipStart_EOM = find(tAnalog >= min(trialArt_ts) + snipRange_EOM(1), 1, 'first');
           snipEnd_EOM = snipStart_EOM + (snipRange_EOM(2)-snipRange_EOM(1)) * adfreq ;
           EOM_snip = EOM_data(:, snipStart_EOM : snipEnd_EOM);
           all_EOM_snips = cat(3, all_EOM_snips, EOM_snip);
           
       end

    end

    snipsEMG{ii} = all_EMG_snips;
    snipsEOM{ii} = all_EOM_snips;
    nArt(ii) = size(all_EMG_snips, 3);
    nTrials(ii) = size(all_EOM_snips, 3);
    

%     figure
%     plot(mean(all_artWaves, 1));

%     yline(2*std(stimTA) + mean(stimTA), 'r');
%     yline(-(2*std(stimTA)) + mean(stimTA), 'r');
    
end

save([save_path 'EMG_EOM_snipInfo2_' date '_' num2str(adfreq) 'Hz.mat'], 'snipsEMG', 'snipsEOM', 'nArt', 'allTargets', 'snipRange_EMG', 'snipRange_EOM', 'adfreq', '-v7.3')

%% Eye Movement Calibration

EOM_data_path = ['F:\RawData\COT_ICMS\' monk '\' cortical_region '\Ripple\'];
EOM_save_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\stimTA\'];

%EOM_date = date;
EOM_date = '20210712';

EOM_ns2_file = [EOM_data_path monk(1) '_' EOM_date '_eyemvmt.ns2'];
EOM_nev_file = [EOM_data_path monk(1) '_' EOM_date '_eyemvmt.nev'];

[marker_ts, marker_data] = FMLReadMarkers(EOM_nev_file);
ts_EOMx = FMLReadTimeseries(EOM_ns2_file, 'analog 5');
ts_EOMy = FMLReadTimeseries(EOM_ns2_file, 'analog 6');

targets = [0, 3, 6, 9];
eyePos = zeros(length(targets),2);

for ti = 1:length(targets)
    
    startTrial = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.');
%     successMarker = Marker(119);
    enterMarker = Marker(248);
    targetMarker = Marker(targets(ti));
    targetPattern = Pattern(startTrial, targetMarker, enterMarker);%, successMarker);
    
    [T_x, ~, PI_x] = TimeseriesFindExtractAlign(ts_EOMx, marker_ts, marker_data, targetPattern,'AlignAtIndex',1,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);
    [T_y, ~, PI_y] = TimeseriesFindExtractAlign(ts_EOMy, marker_ts, marker_data, targetPattern,'AlignAtIndex',1,'ExtractMode','around','ExtractAroundIndex',3,'ExtractBeforeDuration',0.001, 'ExtractAfterDuration', 0.001);

    xCollect = [];
    yCollect = [];
    
    for utsi = 1:length(T_x{1})
        
        xCollect = [xCollect; T_x{1}{utsi}.Data'];
        yCollect = [yCollect; T_y{1}{utsi}.Data'];
        
    end
    
    eyePos(ti, 1) = mean(xCollect);
    eyePos(ti, 2) = mean(yCollect); 
    
end

annuliDist = 2.875;
screenDist = 36.25;

targetVolt = mean(abs(eyePos), 'all');
degPerVolt = atan(annuliDist/screenDist)*360/(2*pi()) / targetVolt;

plot(eyePos(:,1)*degPerVolt, eyePos(:,2)*degPerVolt, 'o');
axis([-10 10 -10 10])

 %% Plotting

load([save_path 'EMG_EOM_snipInfo2_' date '_' num2str(adfreq) 'Hz.mat']);
pltRange_EOM = [0 200]/1000;

EMG_fig = figure('Position', [500 100 900 875]);
EOM_fig = figure('Position', [500 100 900 875]);
EOM_fig2 = figure('Position', [500 100 700 875]);

muscleNames = ["Triceps", "Biceps", "Forearm Flex 1", "Forearm Flex 2", "Forearm Ext 1",...
                "Forearm Ext 2"];%, "Deltoid", "Triceps2"];      
eyeNames = ["EOM X", "EOM Y"];

plotTime_EMG = [snipRange_EMG(1) : 1/adfreq : snipRange_EMG(2)]*1000;
plotTime_EOM = [snipRange_EOM(1) : 1/adfreq : snipRange_EOM(2)]*1000;
    

for ii = 1:length(nArt)
    
    figNums_EMG = [ii:length(allTargets):length(muscleNames)*length(allTargets)];
    figNums_EOM = [ii:length(allTargets):length(eyeNames)*length(allTargets)];
    stimTA_EMG = mean(snipsEMG{ii}, 3);
    stimTrig_EOM = (snipsEOM{ii});
    stimTrig_EOM = stimTrig_EOM(:, 1:length(pltRange_EOM(1):1/adfreq:pltRange_EOM(2)), :);
    plotTime_EOM = plotTime_EOM(:, 1:length(pltRange_EOM(1):1/adfreq:pltRange_EOM(2)), :);
    
    figure(EMG_fig)
    for pi = 1:length(muscleNames)
        
        this_stimTA_EMG = stimTA_EMG(pi,:);
        linReg = fitlm(plotTime_EMG, this_stimTA_EMG);
        if linReg.Coefficients.pValue(2) < 0.05
            flat_stimTA_EMG = this_stimTA_EMG - (plotTime_EMG*linReg.Coefficients.Estimate(2));
        else
            flat_stimTA_EMG = this_stimTA_EMG;
        end

        SD2 = 2*std(flat_stimTA_EMG(plotTime_EMG < 0));
        EMG_Axes{pi, ii} = subplot(length(muscleNames),length(allTargets),figNums_EMG(pi));
        hold on
        EMG_Plot{pi, ii} = plot(plotTime_EMG, flat_stimTA_EMG);
        xline(0, '--k', 'LineWidth', 1.5);
        yline(SD2 + mean(flat_stimTA_EMG), '--r');
        yline(-SD2 + mean(flat_stimTA_EMG), '--r');
        hold off
        xlim(snipRange_EMG*1000)
        if ii == 1
            ylabel(muscleNames(pi))
        end
        if figNums_EMG(pi) == max(figNums_EMG)
            xlabel('Time (msec)')
        end
        if figNums_EMG(pi) == min(figNums_EMG)
            title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(nArt(ii)) ')'])
        end

    end

    figure(EOM_fig)
    for mi = 1:length(eyeNames)

        EOM_Axes{mi, ii} = subplot(length(eyeNames),length(allTargets),figNums_EOM(mi));
        hold on
        for pi2 = 1:size(stimTrig_EOM, 3)        

            EOM_Plot{pi2, ii} = plot(plotTime_EOM, stimTrig_EOM(mi,:,pi2)*degPerVolt, 'color', [0,0,0]+0.8);
%           xline(0, '--k', 'LineWidth', 1.5);
        end
        meanEOM = mean(stimTrig_EOM, 3);
        plot(plotTime_EOM, meanEOM(mi,:)*degPerVolt, 'color', 'k', 'LineWidth', 2.0);
        hold off

%         xlim(snipRange_EOM*1000)
%         xlim(pltRange_EOM*1000)
        if ii == 1
            ylabel(eyeNames(mi))
        end
        if figNums_EOM(mi) == max(figNums_EOM)
            xlabel('Time (msec)')
        end
        if figNums_EOM(mi) == min(figNums_EOM)
            title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(size(snipsEOM{ii}, 3)) ')'])
        end

    end

    figure(EOM_fig2)
    EOM_Axes2{ii} = subplot(length(nArt),1,ii);
    hold on  
    for pi3 = 1:size(stimTrig_EOM, 3)        

        EOM_Plot2{pi3, ii} = plot(stimTrig_EOM(1,:,pi3)*degPerVolt, stimTrig_EOM(2,:,pi3)*degPerVolt, 'color', [0,0,0]+0.9);
    %           xline(0, '--k', 'LineWidth', 1.5);
    end
    meanEOM = mean(stimTrig_EOM, 3);
    plot(meanEOM(1,:)*degPerVolt, meanEOM(2,:)*degPerVolt, 'color', 'k', 'LineWidth', 2.0);
    hold off
%     xlim([-5 5])
%     ylim([-5 5])
    if ii == length(allTargets)
        xlabel('X Position')
    end
    ylabel('Y Position')
    title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(size(snipsEOM{ii}, 3)) ')'])
    
end

figure(EMG_fig)
suptitle('Stimulus-Triggered EMG Averages');
% figure(EOM_fig)
% suptitle('Stimulus-Triggered Eye Movement Averages');
% figure(EOM_fig2)
% suptitle(['Stimulus-Triggered Eye Movements Over ' num2str(pltRange_EOM(2)*1000) ' msec']);

saveas(EMG_fig, [save_path monk(1) '_' date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.png'])
saveas(EMG_fig, [save_path monk(1) '_' date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.fig'])
saveas(EOM_fig, [save_path monk(1) '_' date '_EOMstimTA_Time' M1_elec '_' num2str(adfreq) 'Hz.png'])
saveas(EOM_fig, [save_path monk(1) '_' date '_EOMstimTA_Time' M1_elec '_' num2str(adfreq) 'Hz.fig']) 
saveas(EOM_fig2, [save_path monk(1) '_' date '_EOMstimTA_Pos' M1_elec '_' num2str(adfreq) 'Hz.png'])
saveas(EOM_fig2, [save_path monk(1) '_' date '_EOMstimTA_Pos' M1_elec '_' num2str(adfreq) 'Hz.fig']) 
 