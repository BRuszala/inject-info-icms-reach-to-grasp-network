close all; clear;

%% Initial Settings
addpath(genpath('D:\MATLAB\FML_matlab'))
addpath(genpath('D:\MATLAB\fileHandling_BR'))

monk = 'Felix';
cortical_region = 'PMv';

EMG_date = '20221206';
emgFileEnd = 'EMGbaseline';    %After COT_ICMS
artFileEnd = 'EMG_M1catch_78_art_207_WLsort';    %After COT_ICMS

useNs5 = true;
useNs2 = false;
usePlx = false;
recArtTs = true;
emgToken = [monk(1) '_' EMG_date '_COT_ICMS_' emgFileEnd];
artToken = [monk(1) '_' EMG_date '_COT_ICMS_' artFileEnd];

% task_targets = [0, 2, 4, 6];
task_targets = [0, 1, 2, 3];
EMG_target = [];
% EMG_target = [NaN]; %Use NaN for catch trials
M1_elec = ['78'];
artChan = 207;

snipRange_EMG = [-20 40]/1000; % In sec

% data_path = ['D:\RawData\COT_ICMS\' monk '\' cortical_region '\Ripple\'];
data_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\stimTA\emgData\'];
save_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\stimTA\'];


%% Data Extraction

if ~usePlx
    try
        nev_filename = [data_path emgToken '.nev'];
        [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
    catch
        nev_filename = [data_path emgToken '_shrink.nev'];
        [marker_ts, marker_data] = FMLReadMarkers(nev_filename);
    end
else
    [marker_ts, marker_data] = FMLReadMarkers([data_path emgToken '.plx']);
end

if useNs2
    EMG_fullFile = [data_path emgToken '.ns2'];   
elseif useNs5
    EMG_fullFile = [data_path emgToken '.ns5'];  
elseif usePlx
    EMG_fullFile = [data_path emgToken '.plx'];  
end

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
else
    if ~usePlx
        art_filename = [data_path artToken '.nex'];
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
    else
        
        art_filename = [data_path artToken '.plx'];
        %Reading Plexon file (get markers and Xdsip and Ydisp Time series data)
        [~, ~, artTs_1, wave_1] = plx_waves_v(art_filename, artChan, 1);
        [~, ~, artTs_2, wave_2] = plx_waves_v(art_filename, artChan, 2);
        [~, ~, artTs_3, wave_3] = plx_waves_v(art_filename, artChan, 3);
        [~, ~, artTs_4, wave_4] = plx_waves_v(art_filename, artChan, 4);
            
        artTs = [artTs_1; artTs_2; artTs_3; artTs_4];
    end
end


if ~usePlx
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 13');
    adfreq = 1/tempEMG.Period;
    ts = tempEMG.StartTime;

    %  60 Hz notch filter
    nFilt60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'DesignMethod','butter','SampleRate', adfreq);
        
    EMG_data(1,:) = filtfilt(nFilt60, tempEMG.Data); 
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 14');
    EMG_data(2,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 15');
    EMG_data(3,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 16');
    EMG_data(4,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 17');
    EMG_data(5,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 18');
    EMG_data(6,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 19');
    EMG_data(7,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'analog 20');
    EMG_data(8,:) = filtfilt(nFilt60, tempEMG.Data);
else
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD09');
    adfreq = 1/tempEMG.Period;
    ts = tempEMG.StartTime;

    %  60 Hz notch filter
    nFilt60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,'DesignMethod','butter','SampleRate', adfreq);
    
    EMG_data(1,:) = filtfilt(nFilt60, tempEMG.Data); 
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD10');
    EMG_data(2,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD11');
    EMG_data(3,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD12');
    EMG_data(4,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD13');
    EMG_data(5,:) = filtfilt(nFilt60, tempEMG.Data);
    tempEMG = FMLReadTimeseries(EMG_fullFile, 'AD14');
    EMG_data(6,:) = filtfilt(nFilt60, tempEMG.Data);

end

% rectAD = abs(ad);
EMG_data = abs(EMG_data);


%% Triggered Averaging

figure()

tAnalog = [1/adfreq:1/adfreq:length(EMG_data)/adfreq]' + ts;
allTargets = [task_targets, EMG_target];
numBootstrap = 100;
jitterSigma = 10; %in msec

startMarker = Marker(251, '.', '.', '.', '.', '.', '.', '.', '.', '.');
instructMarker = Marker(210);
endMarker = Marker(253); 
stimDur = cell(length(allTargets),1);
disp('Starting stimTA Calculation...')
for ii = 1:length(allTargets)

    if ~isnan(allTargets(ii))
        targetMarker = Marker(allTargets(ii));
        targetPattern = Pattern(startMarker.replace_at(7, 0), targetMarker, instructMarker, endMarker);
    else
        targetPattern = Pattern(startMarker.replace_at(7, 10), instructMarker, endMarker);
    end

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
    disp(['StimTAs ' num2str(ii/length(allTargets)*100, '%d') '% Done'])
%     figure
%     plot(mean(all_artWaves, 1));

%     yline(2*std(stimTA) + mean(stimTA), 'r');
%     yline(-(2*std(stimTA)) + mean(stimTA), 'r');
    
    
end

% jitterEMG_snip = cell(numBootstrap,length(allTargets);
% for ji = 1:numBootstrap
% 
%     jitterArtTs = artTs + normrnd(0, jitterSigma, 1, length(artTs))/1000;
%     for ii = 1:length(allTargets)
%     
%         if ~isnan(allTargets(ii))
%             targetMarker = Marker(allTargets(ii));
%             targetPattern = Pattern(startMarker.replace_at(7, 0), targetMarker, instructMarker, endMarker);
%         else
%             targetPattern = Pattern(startMarker.replace_at(7, 10), instructMarker, endMarker);
%         end
%     
%         patNDX = FindPattern(marker_data, targetPattern);
%     
%         all_jitter_snips = [];
%         all_artWaves = [];
%         stimDur{ii} = zeros(size(patNDX,1),1);
%         for ti = 1:size(patNDX, 1)
%     
%             trialTime = [marker_ts(patNDX(ti,1)) marker_ts(patNDX(ti,end))];
%             art_ts_ndx = (artTs >= trialTime(1) & artTs <= trialTime(2));
%             trialArt_ts = jitterTs(art_ts_ndx);
%     
%                 for ai = 1:length(trialArt_ts)
%                 
%                     jitterSnipStart = find(tAnalog >= trialArt_ts(ai)+snipRange_EMG(1), 1, 'first');
%                     jitterSnipEnd = jitterSnipStart + (snipRange_EMG(2)-snipRange_EMG(1)) * adfreq;
%                     jitterSnip = EMG_data(:, jitterStart_EMG : jitterEnd_EMG);        
%                     all_jitter_snips = cat(3, all_jitter_snips, jitterSnip);
%                 
%                 end
% 
%         end
%     
%         jitterEMG_snip{ji,ii} = all_jitter_snips;
%         nArt(ii) = size(all_jitter_snips, 3);
% 
%     end
% end

save([save_path 'EMG_snipInfo_' EMG_date '_' num2str(adfreq) 'Hz.mat'], 'snipsEMG', 'nArt', 'allTargets', 'snipRange_EMG', 'adfreq', '-v7.3')

 %% Plotting
close all

monk = 'Qulio';
cortical_region = 'PMv';
EMG_date = '20210414'; 
adfreq = 5000;
M1_elec = ['80_78_77_76'];
yLim = [-7 7];
artWin = [-2 4];
excludeMuscles = logical(1);
muscleNums = [3 6];

snipRange_EMG = [-20 40]/1000; % In sec
smooth = true;
smoothWin = 10;

% targetNames = ["T1", "T2", "T3", "T4", "M1 Catch"];
% targetNames = ["T1", "T2", "T3", "M1 Catch"];
targetNames = ["T1", "T2", "T3", "T4"];

% muscleNames = ["Triceps", "Biceps", "Flexors", "Flexors", "Extensors",...
%                "Extensors", "Triceps", "Deltoid"]; 

muscleNames = ["Triceps", "Biceps", "Flexors", "Flexors", "Extensors",...
               "Extensors"]; 

useMuscle = true(size(muscleNames));
if excludeMuscles
    useMuscle(muscleNums) = false;
end

data_path = ['F:\RawData\COT_ICMS\' monk '\' cortical_region '\Ripple\'];
save_path = ['F:\Projects\COT_ICMS\' monk '\' cortical_region '\StimTA\'];


load([save_path 'EMG_snipInfo_' EMG_date '_' num2str(adfreq) 'Hz.mat']);
EMG_fig = figure('Position', [500 100 700 900]);

num_M1_trig = round(mean(nArt(1:4)));
% if num_M1_trig < nArt(5) %To make M1 triggers and others similar in count
%     snipsEMG{5} = snipsEMG{5}(:, :, randperm(nArt(5), num_M1_trig));
% end

plotTime_EMG = [snipRange_EMG(1) : 1/adfreq : snipRange_EMG(2)]*1000; 
if smooth
    plotTime_EMG = plotTime_EMG(ceil(smoothWin/2):end-ceil(smoothWin/2));
end

zScore_stimTA_EMG = cell(sum(useMuscle), length(nArt));
for ii = 1:length(targetNames)
    
    figNums_EMG = [ii:length(allTargets):length(muscleNames)*length(allTargets)];
    avg_stimTA_EMG = mean(snipsEMG{ii}, 3);

    figure(EMG_fig)
    pi = 1;
    for mi = 1:length(muscleNames)
        
        if useMuscle(mi)
            
            if smooth
                this_stimTA_EMG = conv(avg_stimTA_EMG(mi,:), (1/smoothWin)*ones(1,smoothWin), 'same'); %Filtered
                this_stimTA_EMG = this_stimTA_EMG(ceil(smoothWin/2):end-ceil(smoothWin/2));
            else
                this_stimTA_EMG = avg_stimTA_EMG(mi,:);
            end
            linReg = fitlm(plotTime_EMG, this_stimTA_EMG);
            if linReg.Coefficients.pValue(2) < 0.05
                flat_stimTA_EMG = this_stimTA_EMG - (plotTime_EMG*linReg.Coefficients.Estimate(2));
            else
                flat_stimTA_EMG = this_stimTA_EMG;
            end
            zScore_stimTA_EMG{mi,ii} = (flat_stimTA_EMG - mean(flat_stimTA_EMG(plotTime_EMG < 0))) ./ std(flat_stimTA_EMG(plotTime_EMG < 0));
            
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
                xlabel('Time (ms)')
            end
            if figNums_EMG(pi) == min(figNums_EMG)
                title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(size(snipsEMG{ii}, 3)) ')'])
            end
            pi = pi+1;
            
        end
            
    end

end

%%

EMG_fig2 = figure('Position', [500 100 750 900]);
figure(EMG_fig2)
pltEMG = zScore_stimTA_EMG(useMuscle,:);
pltNames = muscleNames(useMuscle);
sigSTD = 3;

pltOrder = [1 2 3 4];
% pltOrder = [5 4 1 3 2];
% pltOrder = [1, 2, 4, 3]
% pltOrder = [1, 4, 3, 2];
% nArt(5) = num_M1_trig;
for pi = 1:size(pltEMG,1)
    
    EMG_Axes2{pi} = subplot(size(pltEMG,1), 1, pi, 'FontSize', 14);
    hold on
    for ii = 1:length(targetNames)
        
%             EMG_Plot{pi, ii} = plot(plotTime_EMG, flat_stimTA_EMG);
        pltEMG{pltOrder(pi),ii}((plotTime_EMG > artWin(1)) & (plotTime_EMG < artWin(2))) = NaN;
        if any(abs(pltEMG{pltOrder(pi),ii}) > sigSTD)
            aboveSTD = find(diff(pltEMG{pltOrder(pi),ii} > sigSTD)==-1);
            belowSTD = find(diff(pltEMG{pltOrder(pi),ii} > sigSTD)==1);            
            if any(((aboveSTD-belowSTD)/adfreq) > 0.001)
                sigStr = '*';  

            else
                sigStr = '';
            end
        else
            sigStr = '';
        end

        plot(plotTime_EMG, pltEMG{pltOrder(pi),ii}, 'LineWidth', 1.5, 'DisplayName', [char(targetNames(ii)) ' (n = ' num2str(size(snipsEMG{ii}, 3)) ')' sigStr]);
%         if ii == 4
%             plot(plotTime_EMG, pltEMG{pltOrder(pi),ii}, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'DisplayName', [char(targetNames(ii)) ' (n = ' num2str(size(snipsEMG{ii}, 3)) ')' sigStr]);
%         else
%             plot(plotTime_EMG, pltEMG{pltOrder(pi),ii}, 'LineWidth', 2, 'DisplayName', [char(targetNames(ii)) ' (n = ' num2str(size(snipsEMG{ii}, 3)) ')' sigStr]);
%         end
    end
%     xline(0, '-r', 'LineWidth', 1.5);
%             yline(SD2 + mean(flat_stimTA_EMG), '--r');
%             yline(-SD2 + mean(flat_stimTA_EMG), '--r');
    yline(sigSTD , '-.k', 'LineWidth', 1.5);
    yline(-sigSTD , '-.k', 'LineWidth', 1.5);
    yline(0, '-k', 'LineWidth', 2);
    hold off
%             xlim([snipRange_EMG(1)+0.001 snipRange_EMG(2)-0.001]*1000)
    xlim([snipRange_EMG]*1000)
%     xlim([-20 40])
    ylim(yLim)
    patch([artWin flip(artWin)], [yLim(1) yLim(1) yLim(2) yLim(2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    ax = gca;
    box on
        
%     legend(flip(EMG_Axes2{pi}.Children(end-3:end)), 'location', 'eastoutside', 'FontSize', 12)
    legend(flip(EMG_Axes2{pi}.Children(end-4:end)), 'location', 'eastoutside', 'FontSize', 12)
    ylabel(pltNames(pltOrder(pi)))
    if pi == size(pltEMG,1)
        xlabel('Time (ms)')
    end
    ax.Position(3) = 0.5;
    ax.LineWidth = 1.25;
    
%     if pi == 1
%         title(['Target ' num2str(allTargets(ii)) ' (n= ' num2str(nArt(ii)) ')'])
%     end

end

% tb_x = 0.13;
% tb_y = 0.95;
% tb_w = 0.77;
% tb_h = 0.03;
% annStr = ['T1 (n = ' num2str(nArt(1)) '), T2 (n = ' num2str(nArt(2)) '), T3 (n = ' num2str(nArt(3)) '), T4 (n = ' num2str(nArt(4)) '), M1 catch (n = ' num2str(nArt(5)) ')'];  
% annotation('Textbox', [tb_x tb_y tb_w tb_h], 'String', annStr,'EdgeColor','k', 'FontSize', 11)

% legend(flip(EMG_Axes{end}.Children(end-4:end)), 'location', 'southoutside', 'FontSize', 12)

% suptitle('Stimulus-Triggered EMG Averages');
EMG_Axes2{1}.Position = [0.08, 0.82, 0.6, 0.14];
EMG_Axes2{2}.Position = [0.08, 0.63, 0.6, 0.14];
EMG_Axes2{3}.Position = [0.08, 0.45, 0.6, 0.14];
EMG_Axes2{4}.Position = [0.08, 0.26, 0.6, 0.14];
% EMG_Axes2{5}.Position = [0.08, 0.07, 0.6, 0.14];

EMG_Axes2{1}.Legend.Position = [0.7102 0.8283 0.2700 0.1233];
EMG_Axes2{2}.Legend.Position = [0.7102 0.6383 0.2700 0.1233];
EMG_Axes2{3}.Legend.Position = [0.7102 0.4583 0.2700 0.1233];
EMG_Axes2{4}.Legend.Position = [0.7102 0.2683 0.2700 0.1233];
% EMG_Axes2{5}.Legend.Position = [0.7102 0.0783 0.2700 0.1233];

saveas(EMG_fig2, [save_path monk(1) '_' EMG_date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.png'])
% saveas(EMG_fig, [save_path monk(1) '_' EMG_date '_EMGstimTA' M1_elec '_' num2str(adfreq) 'Hz.fig'])


 