clear

clc

[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*.histperFork', 'MatFile data:','MultiSelect', 'on');
coLocData = [];
binning = [0:.1:5];
for i=1:length(ROIDataFilename)
    
    if i==1
        coLocData = importdata(ROIDataFilename{i});
    else
        coLocDatatemp = importdata(ROIDataFilename{i});
        coLocData.minDist = [coLocData.minDist; coLocDatatemp.minDist];
        coLocData.minDistImmobilised = [coLocData.minDistImmobilised; coLocDatatemp.minDistImmobilised];
        coLocData.minDistLoc = [coLocData.minDistLoc; coLocDatatemp.minDistLoc];
        coLocData.minDistImmoLoc = [coLocData.minDistImmoLoc; coLocDatatemp.minDistImmoLoc];
    end
    
    
end

set(0,'DefaultLineMarkerSize',20,'DefaultLineLineWidth',2);
figure,
histC.minDist = histc(coLocData.minDist, binning);
hminDist = plot(binning,histC.minDist/sum(histC.minDist),'-b.');


hold on
histC.minDistImmobilised = histc(coLocData.minDistImmobilised, binning);
hminDistImmobilised = plot(binning,histC.minDistImmobilised/sum(histC.minDistImmobilised),'-r.');

legend([hminDist hminDistImmobilised],...
    ['Sample size = ',num2str(length(coLocData.minDist))],...
    ['Sample size = ',num2str(length(coLocData.minDistImmobilised))])
title('Mean track position distance to DNAQ')

figure,
histC.minDistLoc = histc(coLocData.minDistLoc, binning);
hminDistLoc = plot(binning,histC.minDistLoc/sum(histC.minDistLoc),'-b.');

hold on
histC.minDistImmoLoc = histc(coLocData.minDistImmoLoc, binning);
hminDistImmoLoc = plot(binning,histC.minDistImmoLoc/sum(histC.minDistImmoLoc),'-r.');


legend([hminDistLoc hminDistImmoLoc],...
    ['Sample size = ',num2str(length(coLocData.minDistLoc))],...
    ['Sample size = ',num2str(length(coLocData.minDistImmoLoc))])
title('Track pos distance to DNAQ')

figure,
plot(binning,cumsum(histC.minDist/sum(histC.minDist)),'-b.')
hold on
plot(binning,cumsum(histC.minDistImmobilised/sum(histC.minDistImmobilised)),'-r.')

figure,
plot(binning,cumsum(histC.minDistLoc/sum(histC.minDistLoc)),'-b.')
hold on
plot(binning,cumsum(histC.minDistImmoLoc/sum(histC.minDistImmoLoc)),'-r.')








