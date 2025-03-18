close all; clear;

%% Initial Settings
addpath(genpath('C:\Users\bruszala\Documents\MATLAB\FML_matlab'))
monk = 'Qulio';
cortical_region = 'AIP';

data_path = ['D:\COT_ICMS\' monk '\' cortical_region '\Ripple\'];
save_path = ['C:\Users\bruszala\Documents\Projects\COT_ICMS\' monk '\' cortical_region '\StimTA\'];

EMG_date = '20211210';

useNs5 = true;
useNs2 = false;
rippleToken = ['Q_20211210_COT_ICMS_M1catch_78'];

task_targets = [0, 2, 4, 6];
EMG_target = [1];
M1_elec = '78';

snipRange_EMG = [-30 50]/1000; % In sec

%% Data Extraction

nev_filename = [data_path rippleToken '.nev'];
[marker_ts, marker_data] = FMLReadMarkers(nev_filename);

if useNs2
     EMG_fullFile = [data_path rippleToken '.ns2'];   
elseif useNs5
     EMG_fullFile = [data_path rippleToken '.ns5'];   
end

sma_chan = 2;
artTs = FMLReadNEVDIOEventsTimestamps(nev_filename, sma_chan);
ISI = diff(artTs);
repeatNDX = [false ISI<0.001];
artTs(repeatNDX) = [];

tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 13');
EMG_data(1,:) = tempEMG.Data; 
adfreq = 1/tempEMG.Period;
ts = tempEMG.StartTime;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 14');
EMG_data(2,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 15');
EMG_data(3,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 16');
EMG_data(4,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 17');
EMG_data(5,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 18');
EMG_data(6,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 19');
EMG_data(7,:) = tempEMG.Data;
tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 20');
EMG_data(8,:) = tempEMG.Data;

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
stimDur = cell(length(allTargets),1);

for ii = 1:length(allTargets)

    targetMarker = Marker(allTargets(ii));
    targetPattern = Pattern(startMarker, targetMarker, instructMarker, endMarker);

    patNDX = FindPattern(marker_data, targetPattern);

    all_EMG_snips = [];
    all_artWaves = [];
    stimDur{ii} = zeros(size(patNDX,1),1);
    for ti = 1:size(patNDX, 1)

       trialTime = [marker_ts(patNDX(ti,1)) marker_ts(patNDX(ti,end))];
       art_ts_ndx = (artTs >= trialTime(1) & artTs <= trialTime(2));
       trialArt_ts = artTs(art_ts_ndx);
       
       for ai = 1:length(trialArt_ts)

           snipStart_EMG = find(tAnalog >= trialArt_ts(ai)+snipRange_EMG(1), 1, 'first');
           snipEnd_EMG = snipStart_EMG + (snipRange_EMG(2)-snipRange_EMG(1)) * adfreq;
           
           EMG_snip = EMG_data(:, snipStart_EMG : snipEnd_EMG);
           all_EMG_snips = cat(3, all_EMG_snips, EMG_snip);
           
       end

    end

    snipsEMG{ii} = all_EMG_snips;
    nArt(ii) = size(all_EMG_snips, 3);
   
%     figure
%     plot(mean(all_artWaves, 1));

%     yline(2*std(stimTA) + mean(stimTA), 'r');
%     yline(-(2*std(stimTA)) + mean(stimTA), 'r');
    
end

save([save_path 'EMG_snipInfo_' EMG_date '_' num2str(adfreq) 'Hz.mat'], 'snipsEMG', 'nArt', 'allTargets', 'snipRange_EMG', 'adfreq', '-v7.3')

 %% Plotting
 
adfreq = 30000;
load([save_path 'EMG_snipInfo_' EMG_date '_' num2str(adfreq) 'Hz.mat']);

EMG_fig = figure('Position', [500 100 900 875]);

muscleNames = ["Triceps", "Biceps", "Forearm Flex 1", "Forearm Flex 2", "Forearm Ext 1",...
                "Forearm Ext 2", "Triceps", "Deltoid"]; 
useMuscle = true(size(muscleNames));
useMuscle(1,7) = false;

plotTime_EMG = [snipRange_EMG(1) : 1/adfreq : snipRange_EMG(2)]*1000;  
for ii = 1:length(nArt)
    
    figNums_EMG = [ii:length(allTargets):length(muscleNames)*length(allTargets)];
    stimTA_EMG = mean(snipsEMG{ii}, 3);

    figure(EMG_fig)
    pi = 1;
    for mi = 1:length(muscleNames)
        
        if useMuscle(mi)
            this_stimTA_EMG = stimTA_EMG(mi,:);
            linReg = fitlm(plotTime_EMG, this_stimTA_EMG);
            if linReg.Coefficients.pValue(2) < 0.05
                flat_stimTA_EMG = this_stimTA_EMG - (plotTime_EMG*linReg.Coefficients.Estimate(2));
            else
                flat_stimTA_EMG = this_stimTA_EMG;
            end

            SD2 = 2*std(flat_stimTA_EMG(plotTime_EMG < 0));
            EMG_Axes{pi, ii} = subplot(sum(useMuscle),length(allTargets),figNums_EMG(pi));
            hold on
            EMG_Plot{pi, ii} = plot(plotTime_EMG, flat_stimTA_EMG);
            xline(0, '--k', 'LineWidth', 1.5);
            yline(SD2 + mean(flat_stimTA_EMG), '--r');
            yline(-SD2 + mean(flat_stimTA_EMG), '--r');
            hold off
            xlim(snipRange_EMG*1000)
            if ii == 1
                ylabel(muscleNames(mi))
            end
            if figNums_EMG(pi) == max(figNums_EMG)
                xlabel('Time (msec)')
            end
            if figNums_EMG(pi) == min(figNums_EMG)
                title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(nArt(ii)) ')'])
            end
            pi = pi+1;
        end

    end

end

figure(EMG_fig)
suptitle('Stimulus-Triggered EMG Averages');

saveas(EMG_fig, [save_path monk(1) '_' EMG_date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.png'])
saveas(EMG_fig, [save_path monk(1) '_' EMG_date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.fig'])


 