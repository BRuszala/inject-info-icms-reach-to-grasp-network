clc; clear; close all;

%% Adding paths to FML Toolboxes
addpath(genpath('F:\bruszala\Documents\MATLAB\FML_matlab'))

%% Overlaying all areas for one monkey. Choose File and sorting parameters
close all
monk = 'Felix';
data_path = ['F:\Projects\COT_ICMS\' monk '\'];
savedir = ['F:\Projects\COT_ICMS\'];
corticalRegions = ["S1", "PMv", "AIP"];
useRegions = [true, true, true];
corticalColors = {[250 202 0]/255, [0 150 68]/255, [112 48 160]/255};
paramType = 'durSweep';
sweepType = 'single';    % or 'all' for interleaved

fontSize = 20;
markerSize = 200;
lineWidth = 4;

ciFig = figure();
for ci = 1:length(corticalRegions)
    
    if useRegions(ci)
        ci_sweep = load([data_path '\' char(corticalRegions(ci)) '\ParamSweeps\' monk(1) '_' sweepType 'Fit_' paramType '.mat']); 
        corticalRegions(ci)
        ci_sweep.singleSummary.Model

        if strcmpi(paramType, 'ampSweep')
            titleStr = 'Amplitude Sweep';
            xStr = ['Amplitude (' char(181) 'A)'];
            xRange = [0:0.01:60]';
        elseif strcmpi(paramType, 'freqSweep')
            titleStr = 'Frequency Sweep';
            xStr = ['Frequency (Hz)'];
            xRange = [0:0.01:85]';
        elseif strcmpi(paramType, 'durSweep')
            titleStr = 'Pulse-Train Duration Sweep';
            xStr = ['Pulse-Train Duration (msec)'];
            xRange = [0:0.1:800]';
        end

        if strcmpi(sweepType, 'single')

            if strcmpi(corticalRegions(ci), "S1")
                ci_plotModel = [xRange, sigmoid_fnc(xRange, ci_sweep.singleSummary.Model.a, ci_sweep.singleSummary.Model.b,...
                                                          ci_sweep.singleSummary.Model.c, 33)];
            else
                ci_plotModel = [xRange, sigmoid_fnc(xRange, ci_sweep.singleSummary.Model.a, ci_sweep.singleSummary.Model.b,...
                                                          ci_sweep.singleSummary.Model.c, ci_sweep.singleSummary.Model.d)];
            end
            hold on

            scatter(ci_sweep.singleSummary.Levels, ci_sweep.singleSummary.Performance, markerSize, corticalColors{ci}, 'filled')
            pCI(ci) = plot(ci_plotModel(:,1), ci_plotModel(:,2), 'Color', corticalColors{ci} , 'LineWidth', lineWidth);
            set(gca, 'Fontsize', fontSize)
            pCI(ci).DisplayName = [char(corticalRegions(ci)) ' (R^2 = ' num2str(round(ci_sweep.singleSummary.GoodnessOfFit.rsquare, 3)) ')'];

            hold off

        elseif strcmpi(sweepType, 'all')

            ci_plotModel = [xRange, sigmoid_fnc(xRange, ci_sweep.allSummary.Model.a,...
                                                              ci_sweep.allSummary.Model.b,...
                                                              ci_sweep.allSummary.Model.c, ...
                                                              ci_sweep.singleSummary.Model.d)];
            hold on

            scatter(ci_sweep.allSummary.Levels, ci_sweep.allSummary.Performance, markerSize, corticalColors{ci} , 'filled')
            pCI(ci) = plot(ci_plotModel(:,1), ci_plotModel(:,2), 'Color', corticalColors{ci}, 'LineWidth', lineWidth);
            set(gca, 'Fontsize', fontSize)
            pCI(ci).DisplayName = [char(corticalRegions(ci)) ' (R^2 = ' num2str(round(ci_sweep.allSummary.GoodnessOfFit.rsquare, 3)) ')'];

            hold off

        end
    end
end

legend([pCI(useRegions)], 'location', 'southeast')
ax = gca;
ax.LineWidth = 2;
xlim([xRange(1) xRange(end)])
xlabel(xStr)
ylim([0 100])
ylabel('% Correct')

saveas(gcf, [savedir monk '\Comparisons\' monk(1) '_' paramType 'Comparison_S1_PMv_AIP.png'])

%% Overlaying one area for all monkeys. Choose File and sorting parameters
close all
corticalRegion = 'AIP';
paramType = 'ampSweep';
sweepType = 'single';    % or 'all' for interleaved
savedir = ['F:\Projects\COT_ICMS\'];

fontSize = 18;
markerSize = 80;
lineWidth = 2;

data_path = ['F:\Projects\COT_ICMS\'];
Q_sweep = load([data_path 'Qulio\' corticalRegion '\ParamSweeps\Q_' sweepType 'Fit_' paramType '.mat']); 
F_sweep = load([data_path 'Felix\' corticalRegion '\ParamSweeps\F_' sweepType 'Fit_' paramType '.mat']); 

if strcmpi(paramType, 'ampSweep')
    titleStr = 'Amplitude Sweep';
    xStr = ['Amplitude (' char(181) 'A)'];
    xRange = [0:0.01:60]';
    tickStep = 10;
elseif strcmpi(paramType, 'freqSweep')
    titleStr = 'Frequency Sweep';
    xStr = ['Frequency (Hz)'];
    xRange = [0:0.01:85]';
    tickStep = 10;
elseif strcmpi(paramType, 'durSweep')
    titleStr = 'Pulse-Train Duration Sweep';
    xStr = ['Pulse-Train Duration (msec)'];
    xRange = [0:0.1:800]';
    tickStep = 200;
end

if strcmpi(sweepType, 'single')
    Q_plotModel = [xRange, sigmoid_fnc(xRange, Q_sweep.singleSummary.Model.a,...
                                                      Q_sweep.singleSummary.Model.b,...
                                                      Q_sweep.singleSummary.Model.c, ...
                                                      Q_sweep.singleSummary.Model.d)];
    F_plotModel = [xRange, sigmoid_fnc(xRange, F_sweep.singleSummary.Model.a,...
                                                       F_sweep.singleSummary.Model.b,...
                                                       F_sweep.singleSummary.Model.c, ...
                                                       F_sweep.singleSummary.Model.d)];
    
    figure();
    hold on

    scatter(Q_sweep.singleSummary.Levels, Q_sweep.singleSummary.Performance, markerSize, [0 0.4470 0.7410], 'filled')
    pQ = plot(Q_plotModel(:,1), Q_plotModel(:,2), 'Color', [0 0.4470 0.7410] , 'LineWidth', lineWidth);
    set(gca, 'Fontsize', fontSize)
    Q_lgnd = ['Q (R^2 = ' num2str(round(Q_sweep.singleSummary.GoodnessOfFit.rsquare, 3)) ')'];

    scatter(F_sweep.singleSummary.Levels, F_sweep.singleSummary.Performance, markerSize, [0.8500 0.3250 0.0980], 'filled')
    pF = plot(F_plotModel(:,1), F_plotModel(:,2), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', lineWidth);
    set(gca, 'Fontsize', fontSize)
    F_lgnd = ['F (R^2 = ' num2str(round(F_sweep.singleSummary.GoodnessOfFit.rsquare, 3)) ')'];

    hold off

%     legend([pQ pF], Q_lgnd, F_lgnd, 'location', 'southeast')
    xlim([xRange(1) xRange(end)])
    xticks(xRange(1) : tickStep : xRange(end))
    xlabel(xStr)
    ylim([0 100])
    ylabel('% Correct')
%     title([titleStr ' Across Cortical Areas'])

elseif strcmpi(sweepType, 'all')
    Q_plotModel = [xRange, sigmoid_fnc(xRange, Q_sweep.allSummary.Model.a,...
                                                      Q_sweep.allSummary.Model.b,...
                                                      Q_sweep.allSummary.Model.c, ...
                                                      Q_sweep.singleSummary.Model.d)];
    F_plotModel = [xRange, sigmoid_fnc(xRange, F_sweep.allSummary.Model.a,...
                                                       F_sweep.allSummary.Model.b,...
                                                       F_sweep.allSummary.Model.c, ...
                                                       F_sweep.singleSummary.Model.d)];
    figure
    hold on

    scatter(Q_sweep.allSummary.Levels, 100*Q_sweep.allSummary.Performance, markerSize, [0 0.4470 0.7410], 'filled')
    pQ = plot(Q_plotModel(:,1), 100*Q_plotModel(:,2), 'Color', [0 0.4470 0.7410], 'LineWidth', lineWidth);
    set(gca, 'Fontsize', fontSize)
    Q_lgnd = ['Q (R^2 = ' num2str(round(Q_sweep.allSummary.GoodnessOfFit.rsquare, 3)) ')'];

    scatter(F_sweep.allSummary.Levels, 100*F_sweep.allSummary.Performance, markerSize, [0.8500 0.3250 0.0980], 'filled')
    pF = plot(F_plotModel(:,1), 100*F_plotModel(:,2), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', lineWidth);
    set(gca, 'Fontsize', fontSize)
    F_lgnd = ['F (R^2 = ' num2str(round(F_sweep.allSummary.GoodnessOfFit.rsquare, 3)) ')'];

    hold off

    legend([pQ pF], Q_lgnd, F_lgnd, 'location', 'southeast')
    xlim([xRange(1) xRange(end)])
    xlabel(xStr)
    % xlabel(['Frequency (Hz)'])
    ylim([0 100])
    ylabel('Success Percentage')
%     title([titleStr 'Comparison Across Cortical Areas'])
    
end

Q_sweep.singleSummary.Model
F_sweep.singleSummary.Model
saveas(gcf, [savedir 'MonkComparisons_' paramType '_' sweepType '_AIP.png'])


function ySig = sigmoid_fnc(X, A, B, C,D)

ySig = (A-D)./(1 + exp(-B*(X-C))) + D;

end