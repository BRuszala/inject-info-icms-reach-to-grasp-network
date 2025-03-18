addpath(genpath('F:\MATLAB\Statistics\'))

load('F:\Publications\InjectInformationPMCPPC\Results\rtSuccessComp.mat', 'rtSuccessStruct')

thisMonk = 'Felix';
monkNDX = cellfun(@(x) strcmpi(x, thisMonk), {rtSuccessStruct.Monk});

area1 = 'PMv';
area2 = 'PMd';

area1NDX = cellfun(@(x) strcmpi(x, area1), {rtSuccessStruct.Area});
area2NDX = cellfun(@(x) strcmpi(x, area2), {rtSuccessStruct.Area});

rt1 = rtSuccessStruct(monkNDX & area1NDX).rts;
rt2 = rtSuccessStruct(monkNDX & area2NDX).rts;

% 1000*[nanmean(rt1) nanmean(rt2)]
1000*[nanmedian(rt1) prctile(rt1, 25) prctile(rt1, 75)]
1000*[nanmedian(rt2) prctile(rt2, 25) prctile(rt2, 75)]

% [h, p, ci, stats] = ttest2(rt1, rt2)
[p, h, stats] = ranksum(rt1, rt2)


success1 = sum(rtSuccessStruct(monkNDX & area1NDX).successNums, 'all');
trials1 =  sum(rtSuccessStruct(monkNDX & area1NDX).trialNums, 'all');
error1 = trials1-success1;
success2 = sum(rtSuccessStruct(monkNDX & area2NDX).successNums, 'all');
trials2 =  sum(rtSuccessStruct(monkNDX & area2NDX).trialNums, 'all');
error2 = trials2-success2;

T = [success1 success2; error1 error2]
[h, p, stats] = chiSquare_BR(T, 0.05)


%% Global tests
load('F:\Publications\InjectInformationPMCPPC\Results\rtSuccessComp.mat', 'rtSuccessStruct')

thisMonk = 'Felix';
monkNDX = cellfun(@(x) strcmpi(x, thisMonk), {rtSuccessStruct.Monk});

allAreas = ["S1", "PMv", "PMd", "AIP", "dPPC"];

performT = zeros(3,5);

allRts = [];
rtGrps = [];
cellRts = cell(1,5);

faceColors = {[255 196 31]/255, [0 176 80]/255, [253 142 51]/255, [159 95 207]/255, [255 0 0]/255};
xLabels = {'S1', 'PMv', 'PMd', 'AIP', 'dPPC'};

for ai = 1:length(allAreas)
    
    areaNDX = cellfun(@(x) strcmpi(x, allAreas(ai)), {rtSuccessStruct.Area});
    success = sum(rtSuccessStruct(monkNDX & areaNDX).successNums, 'all');
    numTrials = sum(rtSuccessStruct(monkNDX & areaNDX).trialNums, 'all');
    error = numTrials-success;

    performT(1,ai) = success;
    performT(2,ai) = error;
    performT(3,ai) = 100*success/numTrials;

    theseRts = rtSuccessStruct(monkNDX & areaNDX).rts;
    thisGrp = ai*ones(size(theseRts));

    allRts = [allRts; theseRts];
    rtGrps = [rtGrps; thisGrp];
    cellRts{1,ai} = 1000*theseRts;

end

medRts = cellfun(@(x) round(median(x)), cellRts)
iqrRts = cellfun(@(x) round([ prctile(x, 25), prctile(x, 75), min(x), max(x)]), cellRts, 'uni', 0)

figure()
plt = violin(cellRts, xLabels, 'edgecolor', 'none', 'facealpha', 0.35, 'plotlegend', 0);
for clr = 1:5
    plt(clr).FaceColor = faceColors{clr};
end
ylim([-500 4500])
xticks(1:5)
xticklabels(xLabels)
yticks([-500 0:1000:4000])
ax = gca;
ax.LineWidth = 1;
savePath = ['F:\Publications\InjectInformationPMCPPC\Results\'];
print([savePath  thisMonk '_' 'rtcompKeySessions.tiff'], '-dtiff', '-r300')


[h, p, chiT] = chiSquare_BR(performT(1:2, :), 0.05)

[p, tbl, stats] = kruskalwallis(allRts,rtGrps)

c = multcompare(stats, "criticalvaluetype", "dunn-sidak")

for area1 = 1:5
    for area2 = 1:5
        if area1 ~= area2
            area1
            area2
            [~, p, thisT] = chiSquare_BR(chiT(:,[area1, area2]), 0.05)
        end
    end
end



%%

load('F:\Projects\COT_ICMS\catchTrialRts.mat')

for ci = 1:length(catchInfo)

    catchInfo(ci).monk
    catchInfo(ci).corticalRegion

    catchRts = catchInfo(ci).rts;
    nonCatchRts = catchInfo(ci).nonCatchRts;

    catchRts(isnan(catchRts)) = [];
    nonCatchRts(isnan(nonCatchRts)) = [];

    median(catchRts)
    median(nonCatchRts)

    [p, h] = ranksum(catchRts, nonCatchRts)

    % figure()
    % hold on
    % histogram(catchRts, [0:0.01:1])
    % histogram(nonCatchRts, [0:0.01:1])

end




