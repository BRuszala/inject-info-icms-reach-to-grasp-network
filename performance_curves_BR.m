clc; clear; close all;

%% Adding paths to FML Toolboxes
addpath(genpath('D:\MATLAB\FML_matlab'))


%% Choose File and sorting parameters
monk = 'Qulio';
brain_region = 'PMv';

useHardDrive = false;

if useHardDrive
    data_path = "F:\Qulio Backup\FilesToRun\"; %Note input data_path as string here
    % Place any files to use in a "FilesToRun" folder on the external hard
    % drive in 'Ripple', 'Plexon', and 'UTCS' folders
    % Note, if necessary files span multiple hard drives, do the same thing
    % on each separately and input both drives as string into data_path ["Drive1", "Drive2"];
else
    data_path = ['D:\RawData\COT_ICMS\' monk '\' brain_region '\']; 
end

save_path = ['D:\Projects\COT_ICMS\' monk '\' brain_region '\ParamSweeps\'];  %For saving figure

param = 'durSweep';  %Filename after _ICMS_

includeSingle = true;
paramLevels = [50 100 150 200 300 400 500 750];  %Levels that were swept
% annuliLevels = [0 1 2 3 4];  % Annuli number for corresponding levels

includeAll = false;
allLevels = [2 5 10 15 20 25 30];
% blocksAll = [0 0 1 1 2 2 3 3];
annuliAll = [4 1 5 2 6 3 0 7];

usePriorSuccess = true; 
trialsUse = NaN; %[400:10000]; %Choose trials or use NaN for all trials
                     
                   
%% Extract trial information from .h5 and .plx file
    % note that the following scripts were designed for UTCS_version 1.0.6
        % changes in .h5 structure may apply
if includeSingle
    
    for pi = 1:length(paramLevels)
    
        if useHardDrive
            if strcmpi(brain_region, 'S1')
                for dpi = 1:length(data_path)
                    plx_fullFiles = fullfile([char(data_path(dpi)) 'Plexon\'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.plx']);
                    plx_files{dpi} = dir(plx_fullFiles);
                    plx_fileList = cat(1, plx_files{:});
                end
                numFiles = length(plx_fileList);
            else
                for dpi = 1:length(data_path)
                    nev_fullFiles = fullfile([char(data_path(dpi)) 'Ripple\'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.nev']);
                    nev_files{dpi} = dir(nev_fullFiles);
                    nev_fileList = cat(1, nev_files{:});
                    ns2_fullFiles = fullfile([char(data_path(dpi)) 'Ripple\'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.ns2']);
                    ns2_files{dpi} = dir(ns2_fullFiles);
                    ns2_fileList = cat(1, ns2_files{:});
                end
                numFiles = length(nev_fileList);
            end
            
            for dpi = 1:length(data_path)
                h5_fullFiles = fullfile([char(data_path(dpi)) 'UTCS\'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.h5']);
                h5_files{dpi} = dir(h5_fullFiles);
                h5_fileList = cat(1, h5_files{:});
            end
            
        else            
            if strcmpi(brain_region, 'S2') && strcmpi(monk(1), 'Q')
                plx_fullFiles = fullfile([data_path 'Plexon\'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.plx']);
                plx_fileList = dir(plx_fullFiles);
                numFiles = length(plx_fileList);
            else
                nev_fullFiles = fullfile([data_path 'Ripple\ParamSweeps'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '*.nev']);
                nev_fileList = dir(nev_fullFiles);
                fileNDX = cellfun(@(x) contains(x, [num2str(paramLevels(pi)) '.nev'])|contains(x, [num2str(paramLevels(pi)) '_shrink']),  {nev_fileList(:).name}');
                nev_fileList = nev_fileList(fileNDX);

                ns2_fullFiles = fullfile([data_path 'Ripple\ParamSweeps'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.ns2']);
                ns2_fileList = dir(ns2_fullFiles);
                numFiles = length(nev_fileList);
            end
             
            h5_fullFiles = fullfile([data_path 'UTCS\ParamSweeps'], ['*_ICMS_' param '_' num2str(paramLevels(pi)) '.h5']);
            h5_fileList = dir(h5_fullFiles);  
        end

        for fi = 1:numFiles  %If recordings for a specific level span more than 1 day
        
    %         plx_filename = [data_path 'Plexon\' plx_fileList(fi).name];
            h5_filename = [h5_fileList(fi).folder '\' h5_fileList(fi).name];

            %Reading H5 file
            h5_file_info = h5info(h5_filename);
            trial_data = struct();
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

            % Get trial serial codes and annuli/success information from all trials
            %Can get other information here as well if desired
            numTrials = length(trial_data);
            trialAnnuli = zeros(numTrials, 1);
            trialSuccess = zeros(numTrials, 1);
            trialPriorSuccess = zeros(numTrials, 1);
    %         trialSerial = zeros(numTrials, 4);
            trialKeep = zeros(numTrials, 1);
            catchTrial = zeros(numTrials, 1);

            for ti = 1:numTrials
                startCodeNDX = find(trial_data(ti).event_codes == 251, 1, 'first');
    %             trialSerial(ti, :) = trial_data(ti).event_codes(startCodeNDX+1 : startCodeNDX+4);
                trialSuccess(ti) = trial_data(ti).success;
                trialPriorSuccess(ti) = trial_data(ti).prior_success;
                trialKeep(ti) = sum(trialsUse == ti);
                trialAnnuli(ti) = trial_data(ti).selected_annulus_index;
                catchTrial(ti) = trial_data(ti).catch_trial;
            end

            if isnan(trialsUse)
                trialKeep(:) = true;
            end

            trialAnnuli(~trialKeep) = [];
            trialSuccess(~trialKeep) = [];
            trialPriorSuccess(~trialKeep) = [];
            catchTrial(~trialKeep) = [];
    %         trialSerial(~trialKeep) = [];
            if usePriorSuccess
                trialSuccess(~trialPriorSuccess) = [];
                trialAnnuli(~trialPriorSuccess) = [];
                catchTrial(~trialPriorSuccess) = [];
    %         trialSerial(~trialPriorSuccess) = [];
            end

            if sum(catchTrial > 0)
                catchSuccess = trialSucess(catchTrial);
                trialSuccess(catchTrial) = [];          
            end


            %Ensures that only trials with correct annuli are chosen
    %         trialsNum(fi,pi) = sum(trialAnnuli == annuliLevels(pi));
    %         successNum(fi,pi) = sum(trialSuccess(trialAnnuli == annuliLevels(pi)));
            trialsNum(fi,pi) = length(trialSuccess);
            successNum(fi,pi) = sum(trialSuccess);

        end    
    end
end
%% For all parameters
if includeAll
    
    if useHardDrive
        if strcmpi(brain_region, 'S1') && strcmpi(monk(1), 'Q')
            for dpi = 1:length(data_path)
                plx_fullFiles = fullfile([char(data_path(dpi)) 'Plexon\'], ['*_ICMS_' param '_all*.plx']);
                plx_files{dpi} = dir(plx_fullFiles);
                plx_fileList = cat(1, plx_files{:});
            end
            numFiles = length(plx_fileList);
        else
            for dpi = 1:length(data_path)
                nev_fullFiles = fullfile([char(data_path(dpi)) 'Ripple\'], ['*_ICMS_' param '_all*.nev']);
                nev_files{dpi} = dir(nev_fullFiles);
                nev_fileList = cat(1, nev_files{:});
                ns2_fullFiles = fullfile([char(data_path(dpi)) 'Ripple\'], ['*_ICMS_' param '_all*.ns2']);
                ns2_files{dpi} = dir(ns2_fullFiles);
                ns2_fileList = cat(1, ns2_files{:});
            end
            numFiles = length(nev_fileList);
        end

        for dpi = 1:length(data_path)
            h5_fullFiles = fullfile([char(data_path(dpi)) 'UTCS\'], ['*_ICMS_' param '_all*.h5']);
            h5_files{dpi} = dir(h5_fullFiles);
            h5_fileList = cat(1, h5_files{:});
        end

    else
        if strcmpi(brain_region, 'S1') && strcmpi(monk(1), 'Q')
            plx_fullFiles = fullfile([data_path 'Plexon\'], ['*_ICMS_' param '_all*.plx']);
            plx_fileList = dir(plx_fullFiles);
            numFiles = length(plx_fileList);
        else
            nev_fullFiles = fullfile([data_path 'Ripple\'], ['*_ICMS_' param '_all*.nev']);
            nev_fileList = dir(nev_fullFiles);

            ns2_fullFiles = fullfile([data_path 'Ripple\'], ['*_ICMS_' param '_all*.ns2']);
            ns2_fileList = dir(ns2_fullFiles);
            numFiles = length(nev_fileList);
        end

        h5_fullFiles = fullfile([data_path 'UTCS\'], ['*_ICMS_' param '_all*.h5']);
        h5_fileList = dir(h5_fullFiles);        
    end
    
    trialsNum_all = zeros(numFiles,length(allLevels));
    successNum_all = zeros(numFiles,length(allLevels));
    
    for fi = 1:numFiles  %If recordings for a specific level span more than 1 day
        
%         plx_filename = [data_path 'Plexon\' plx_fileList(fi).name];
        h5_filename = [h5_fileList(fi).folder '\' h5_fileList(fi).name];

        %Reading H5 file
        h5_file_info = h5info(h5_filename);
        trial_data = struct();
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

        % Get trial serial codes and annuli/success information from all trials
        %Can get other information here as well if desired
        numTrials = length(trial_data);
        trialAnnuli = zeros(numTrials, 1);
        trialSuccess = zeros(numTrials, 1);
        trialPriorSuccess = zeros(numTrials, 1);
        trialBlocks = zeros(numTrials, 1);
%         trialSerial = zeros(numTrials, 4);
        trialKeep = zeros(numTrials, 1);

        for ti = 1:numTrials
            startCodeNDX = find(trial_data(ti).event_codes == 251, 1, 'first');
%             trialSerial(ti, :) = trial_data(ti).event_codes(startCodeNDX+1 : startCodeNDX+4);
            trialSuccess(ti) = trial_data(ti).success;
            trialPriorSuccess(ti) = trial_data(ti).prior_success;
            trialKeep(ti) = sum(trialsUse == ti);
            trialAnnuli(ti) = trial_data(ti).selected_annulus_index;
            trialBlocks(ti) = trial_data(ti).instruction_block_index;
        end

        if isnan(trialsUse)
            trialKeep(:) = true;
        end

        trialAnnuli(~trialKeep) = [];
        trialBlocks(~trialKeep) = [];
        trialSuccess(~trialKeep) = [];
        trialPriorSuccess(~trialKeep) = [];
%         trialSerial(~trialKeep) = [];
        if usePriorSuccess
            trialSuccess(~trialPriorSuccess) = [];
            trialAnnuli(~trialPriorSuccess) = [];
            trialBlocks(~trialPriorSuccess) = [];
%         trialSerial(~trialPriorSuccess) = [];
        end

        startNDX = strfind(h5_filename, 'all_')+4;
        endNDX = strfind(h5_filename, '.h5')-1;
        cellLevels = regexp(h5_filename(startNDX:endNDX), '\_', 'split');
        theseLevels = cellfun(@(x) str2double(x), cellLevels);
        if strcmpi(param, 'durSweep')                      
            for pi2 = 1:length(theseLevels)  
                levelNDX = allLevels == theseLevels(pi2);
                trialsNum_all(fi,levelNDX) = trialsNum_all(fi,levelNDX) + sum(trialBlocks == blocksAll(levelNDX));
                successNum_all(fi,levelNDX) =  successNum_all(fi,levelNDX) + sum(trialSuccess(trialBlocks == blocksAll(levelNDX)));
            end
        else
             for pi2 = 1:length(theseLevels)
                levelNDX = allLevels == theseLevels(pi2);
                trialsNum_all(fi,levelNDX) = trialsNum_all(fi,levelNDX) + sum(trialAnnuli == annuliAll(levelNDX));
                successNum_all(fi,levelNDX) =  successNum_all(fi,levelNDX) + sum(trialSuccess(trialAnnuli == annuliAll(levelNDX)));
             end
        end
    
    end
    
end


%% Generate performance curve and fit data
font_size = 14;
chance = 33;
initParams = [77 0.2 13 33]; %May need to adjust for better fit

paramLevels = singleSummary.Levels
perform = singleSummary.Performance
param = singleSummary.Parameter


% sigmoidFit = fittype(@(a, b, c, d, x) (a./(1 + exp(-b*(x-c))))+d,...
%               'independent', {'x'}, 'coefficients', {'a', 'b', 'c', 'd'});

sigmoidFit = fittype(@(a, b, c, d, x) d+((a-d)./(1 + exp(-b*(x-c)))),...
              'independent', {'x'}, 'coefficients', {'a', 'b', 'c', 'd'});

fitOptions = fitoptions(sigmoidFit);
fitOptions.Lower = [0, 0, -Inf, chance];
fitOptions.Upper = [100, Inf, Inf, chance];
fitOptions.StartPoint = initParams;

xFit = [0 : 0.01 : max(paramLevels)+10]';
if includeAll
    xFit = [0 : 0.01 : max([max(paramLevels), max(allLevels)])+10]';
end

pltParam = figure(1);

if includeSingle
    
%     perform = 100*sum(successNum,1)./sum(trialsNum,1);
%     perform(5) = 0.835
%     paramLevels(7) = [400]
%     perform(5) = [0.83]
%     paramLevels(5) = [400]
%     perform(6) = [0.3]
%     paramLevels(6) = [50]


    [fitted_curve, gof] = fit(paramLevels',perform',sigmoidFit,fitOptions)
%     [fitted_curve, gof] = fit(paramLevels',perform',sigmoidFit,'StartPoint', initParams)
    % CI1_95 = confint(fitted_curve);

    figure(pltParam)

    hold on
    p1_color = [0 0 1];
    scatter(paramLevels, perform, 50, p1_color, 'filled')
    p1 = plot(xFit, fitted_curve(xFit), 'Color', p1_color,'LineWidth',2);
    % p1_CI_1 = plot(xFit, sigmoid_fnc(xFit, CI1_95(1,1), CI1_95(1,2), CI1_95(1,3)), 'Color', p1_color,'LineWidth',1);
    % p1_CI_2 = plot(xFit, sigmoid_fnc(xFit, CI1_95(2,1), CI1_95(2,2), CI1_95(2,3)), 'Color', p1_color,'LineWidth',1);
    % p1_patch1 = [xFit; sigmoid_fnc(xFit, CI1_95(1,1), CI1_95(1,2), CI1_95(1,3))]';
    % p1_patch2 = [xFit; sigmoid_fnc(xFit, CI1_95(2,1), CI1_95(2,2), CI1_95(2,3))]';
    % p1_patch1 = [xFit; sigmoid_fnc(xFit, CI1_95(1,1), fitted_curve.b, fitted_curve.c)]';
    % p1_patch2 = [xFit; sigmoid_fnc(xFit, CI1_95(2,1), fitted_curve.b, fitted_curve.c)]';
    % patch('Faces', [1:2*length(xFit)], 'Vertices', [p1_patch1; flip(p1_patch2)], 'FaceColor', p1_color, 'FaceAlpha', 0.2, 'Edgecolor', 'none')
    gof_txt = ['(R^2 = ' num2str(round(gof.rsquare, 3)) ')'];
    legend(p1, ['Single Level per Session ' gof_txt], 'location', 'southeast', 'FontSize', font_size);
    hold off
    ylim([0 100])
    pltParam.CurrentAxes.FontSize = font_size;
    figstr = '';

    singleSummary.Region = brain_region;
    singleSummary.Parameter = param;
    singleSummary.Levels = paramLevels;
    singleSummary.Performance = perform;
    singleSummary.Model = fitted_curve;
    singleSummary.GoodnessOfFit = gof;
    singleSummary.plotModel = [xFit, fitted_curve(xFit)];
    save([save_path monk(1) '_singleFit_' param '.mat'], 'singleSummary')
    
end

if includeAll
    
    perform_all = sum(successNum_all,1)./sum(trialsNum_all,1);
   
%     [fitted_curve_all, gof_all] = fit(allLevels',perform_all',sigmoidFit,'StartPoint',initParams)
    [fitted_curve_all, gof_all] = fit(allLevels',perform_all',sigmoidFit,fitOptions)
    CI2_95 = confint(fitted_curve_all);
    
    figure(pltParam)
    hold on
    p2_color = [0.8500 0.3250 0.0980];
    scatter(allLevels, perform_all, 50, p2_color, 'filled')
    p2 = plot(xFit, fitted_curve_all(xFit), 'Color', p2_color,'LineWidth', 2);
%     p2_CI_1 = plot(xFit, sigmoid_fnc(xFit, CI2_95(1,1), CI2_95(1,2), CI2_95(1,3)), 'Color', p2_color,'LineWidth',1);
%     p2_CI_2 = plot(xFit, sigmoid_fnc(xFit, CI2_95(2,1), CI2_95(2,2), CI2_95(2,3)), 'Color', p2_color,'LineWidth',1);    
%     p2_patch1 = [xFit; sigmoid_fnc(xFit, CI2_95(1,1), CI2_95(1,2), CI2_95(1,3))]';
%     p2_patch2 = [xFit; sigmoid_fnc(xFit, CI2_95(2,1), CI2_95(2,2), CI2_95(2,3))]';
%     p2_patch1 = [xFit; sigmoid_fnc(xFit, CI2_95(1,1), CI2_95(1,2), CI2_95(1,3))]';
%     p2_patch2 = [xFit; sigmoid_fnc(xFit, CI2_95(2,1), CI2_95(2,2), CI2_95(2,3))]';
%     patch('Faces', [1:2*length(xFit)], 'Vertices', [p2_patch1; flip(p2_patch2)], 'FaceColor', p2_color, 'FaceAlpha', 0.2, 'Edgecolor', 'none')   
    gof_all_txt = ['(R^2 = ' num2str(round(gof_all.rsquare, 3)) ')'];
    hold off
    %legend([p1, p2], ['One Level per Session ' gof_txt], ['All Levels Randomized ' gof_all_txt], 'location', 'southeast', 'FontSize', font_size);
    legend(p2, ['All Levels Randomized ' gof_all_txt], 'location', 'southeast', 'FontSize', font_size);
    figstr = '_withAll';
    
    allSummary.Region = brain_region;
    allSummary.Parameter = param;
    allSummary.Levels = allLevels;
    allSummary.Performance = perform_all;
    allSummary.Model = fitted_curve_all;
    allSummary.GoodnessOfFit = gof_all;
    allSummary.plotModel = [xFit, fitted_curve_all(xFit)];
    save([save_path monk(1) '_allFit_' param '.mat'], 'allSummary')
    
end

figure(pltParam)
ylim([0 100])
xlim([min(xFit) max(xFit)])
ylabel('Success Percentage', 'FontSize', font_size)
if includeSingle && includeAll
    legend([p1, p2], ['Single Level per Session ' gof_txt], ['Levels Randomly Interleaved ' gof_all_txt],...
        'location', 'southeast', 'FontSize', font_size-2);
end


if strcmpi(param, 'ampSweep')
%     title('Changing Parameter: ICMS Amplitude')
    xlabel(['Amplitude (' char(181) 'A)'], 'FontSize', font_size)
    saveas(pltParam, [save_path monk(1) '_ICMS_ampSweep' figstr '.png'])
elseif strcmpi(param, 'freqSweep')
%     title('Changing Parameter: ICMS Frequency')
    xlabel('Frequency (Hz)', 'FontSize', font_size)
    saveas(pltParam, [save_path monk(1) '_ICMS_freqSweep' figstr '.png'])
elseif strcmpi(param, 'durSweep')
%     title('Changing Parameter: ICMS Pulse-Train Duration')
    xlabel('Pulse-Train Duration (msec)', 'FontSize', font_size)
    saveas(pltParam, [save_path monk(1) '_ICMS_durSweep' figstr '.png'])
end


% function ySig = sigmoid_fnc(X, A, B, C)
% 
% ySig = A./(1 + exp(-B*(X-C)));
% 
% end
