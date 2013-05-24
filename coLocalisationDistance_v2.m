function coLocalisationDistance_v2( pos, ROIData, ROIDataFilename, params )
%COLOCALISATIONDISTANCE Calculates the colocalised distances of tracks and
%position
%   Calculate the distances of the inserted tracks against the position.

%Version 2 contains MSD1 and MSD2 method of serparating immoblised
%molecules from the mobile. 


set(0,'DefaultLineMarkerSize',5,'DefaultLineLineWidth',2);

cmap = colormap(jet);
pixel = params.pixel;
dT = params.dT;
sigmaNoise = params.sigmaNoise;
threshold = params.DhistThresh;

threshold1 = params.DThresh(1);
threshold2 = params.DThresh(2);
minSteps = params.DhistMinSteps;
rangeD = params.rangeD;
holdFigure = params.holdFigureCheckbox;
minDist = [];
minDistImmobilised = [];
minDistImmoLoc = [];
minDistLoc = [];

% Binning of the histograms
binning = [0:.05:1];

neg = [-1 -1; 1 1];


%D = zeros(nMolecules,1);
%meanPos = zeros(nMolecules,2);

ll = 1;

immobileTracks = [];

%Code from countMoleculesinROIs.m
if holdFigure
    figureHandle = params.figureHandle;
end

if holdFigure
    figure(figureHandle);
    delete(findall(figureHandle,'Type','hggroup','HandleVisibility','off'));
else
    figure;
end
hold all
axis image;

kk = 1;
doubleForkCount = 0;
totalForks = 0;
countDiscarded = 0;

for jj = 1:length(ROIData)
    
    impolyVertices = ROIData(jj,1).impolyVertices;
    posINindex = inpolygon(pos(:,1),pos(:,2),impolyVertices(:,1),impolyVertices(:,2));
    
    if any(posINindex)
        posINtest = pos(posINindex,:);
        %Line test
        
        plot(ROIData(jj,1).imlineVertices(:,1),ROIData(jj,1).imlineVertices(:,2),'c')
        hold on
        
        xProj = 0;
        yProj = 0;
        edgeFork = false;
        
        for posNum=1:size(posINtest,1)
            [xProj(posNum), yProj(posNum)] = projectionToLine(ROIData(jj,1).imlineVertices,posINtest(posNum,:));
            
            lineLength = sqrt(diff(ROIData(jj,1).imlineVertices(:,1))^2+...
                diff(ROIData(jj,1).imlineVertices(:,2))^2);
            
            x0 = min(ROIData(jj,1).imlineVertices(:,1));
            y0 = min(ROIData(jj,1).imlineVertices(:,2));
            
            rProj = sqrt((yProj(posNum)-y0)^2+(xProj(posNum)-x0)^2);
            edgeFork(posNum) = (rProj<2*lineLength/5 || rProj>3*lineLength/5);
            
            if edgeFork(posNum)
                totalForks = totalForks + 1;
                plot(posINtest(posNum,1),posINtest(posNum,2),'cX','MarkerSize',10)
                hold on
            else
                countDiscarded = countDiscarded + 1;
                plot(posINtest(posNum,1),posINtest(posNum,2),'rX','MarkerSize',10)
            end
        end
    end
    
    posIN = pos(posINindex,:)*pixel;
    
    if sum(posINindex)>2
        doubleForkCount=doubleForkCount+1;
    end
    
    tracksIN = ROIData(jj,1).tracks;
    nMolecules = max(tracksIN(:,4));
    tracksIN(:,1:2) = tracksIN(:,1:2)*pixel;
    plot(pos(:,1), pos(:,2),'ro','MarkerSize',10)
    
    if any(posINindex)
        posIN = posIN(edgeFork,:);
    end
    
    if any(~isempty(posIN) && any(edgeFork))
        
        
        for ii = 1:nMolecules
            tempDistImmoLoc = [];
            tempDistLoc = [];
            
            xx = find(tracksIN(:,4)==ii);
            
            if numel(xx)>minSteps
                
                MSD1 = 0; % reset MSD
                MSD2 = 0;
                
                % sum all squared displacement
                for mm = 1:numel(xx)-1
                    
                    MSD1 = MSD1 + ((tracksIN(xx(mm+1),1) - tracksIN(xx(mm),1))^2 +...
                        (tracksIN(xx(mm+1),2) - tracksIN(xx(mm),2))^2);
                    
                end
                
                MSD1 = MSD1/mm; % mean square displacement
                
                for mm = 1:numel(xx)-2
                    
                    MSD2 = MSD2 + ((tracksIN(xx(mm+2),1) - tracksIN(xx(mm),1))^2 +...
                        (tracksIN(xx(mm+2),2) - tracksIN(xx(mm),2))^2);
                    
                end
                
                MSD2 = MSD2/mm; % mean square displacement
                
                % calculate D from MSD and correct for localization noise
                D1(kk) = MSD1/(4*dT) - sigmaNoise^2*pixel^2/dT;
                D2(kk) = MSD2/(8*dT) - sigmaNoise^2*pixel^2/dT;
                
                if D1(kk) < threshold1 && D2(kk) < threshold2
                                       
                    immobileTracks = [immobileTracks; tracksIN(xx,:)];
                    
                    meanPosImmobile(ll,:) = [mean(tracksIN(xx,1)),mean(tracksIN(xx,2))];
                    
                    tempDistImmo = sqrt((posIN(:,1)-meanPosImmobile(ll,1)).^2+(posIN(:,2)-meanPosImmobile(ll,2)).^2);
                    minDistImmobilised = [minDistImmobilised; min(tempDistImmo)];
                    perFork.minDistImmobilised = [minDistImmobilised; tempDistImmo];
                    [~, indexMin] = min(tempDistImmo);
                    
                    plot([posIN(indexMin,1) meanPosImmobile(ll,1)]/pixel,...
                        [posIN(indexMin,2) meanPosImmobile(ll,2)]/pixel,...
                        'r-')
                    plot(meanPosImmobile(ll,1)/pixel,meanPosImmobile(ll,2)/pixel,'r.')
                    
                    %Localisations distances
                    for numYPET=1:size(posIN,1)
                        tempDistImmoLoc(:,numYPET) = sqrt((posIN(numYPET,1)-tracksIN(xx,1)).^2+(posIN(numYPET,2)-tracksIN(xx,2)).^2);
                        perFork.minDistImmoLoc = [minDistImmoLoc; tempDistImmoLoc(:,numYPET)];
                    end
                    minDistImmoLoc = [minDistImmoLoc; min(tempDistImmoLoc,[],2)];
                    
                    cellLength.minDistImmobilised = minDistImmobilised/sqrt(sum(sum(ROIData(jj,1).imlineVertices.*neg*pixel,1).^2));
                    cellLength.minDistImmoLoc = minDistImmoLoc/sqrt(sum(sum(ROIData(jj,1).imlineVertices.*neg*pixel,1).^2));
                    ll = ll + 1;
                    
                else
                    
                    meanPos(ll,:) = [mean(tracksIN(xx,1)),mean(tracksIN(xx,2))];
                    tempDist = sqrt((posIN(:,1)-meanPos(ll,1)).^2+(posIN(:,2)-meanPos(ll,2)).^2);
                    minDist = [minDist; min(tempDist)];
                    perFork.minDist = [minDist; tempDist];
                    [~, indexMin] = min(tempDist);
                    %
                    %                     plot([posIN(indexMin,1) meanPos(ll,1)]/pixel,...
                    %                         [posIN(indexMin,2) meanPos(ll,2)]/pixel,...
                    %                         'b-')
                    %
                    %                     plot(meanPos(ll,1)/pixel,meanPos(ll,2)/pixel,'b.')
                    
                    %Localisations distances
                    for numYPET=1:size(posIN,1)
                        tempDistLoc(:,numYPET) = sqrt((posIN(numYPET,1)-tracksIN(xx,1)).^2+(posIN(numYPET,2)-tracksIN(xx,2)).^2);
                        perFork.minDistLoc = [minDistLoc; tempDistLoc(:,numYPET)];
                    end
                    minDistLoc = [minDistLoc; min(tempDistLoc,[],2)];
                    cellLength.minDist = minDist/sqrt(sum(sum(ROIData(jj,1).imlineVertices.*neg*pixel,1).^2));
                    cellLength.minDistLoc = minDistLoc/sqrt(sum(sum(ROIData(jj,1).imlineVertices.*neg*pixel,1).^2));
                    ll = ll + 1;
                end
                
                
                % reset MSD and update index
                
                kk = kk + 1;
            end
            
        end
    end
    
end


coLocData.minDist = minDist;
coLocData.minDistImmobilised = minDistImmobilised;

coLocData.minDistLoc = minDistLoc;
coLocData.minDistImmoLoc = minDistImmoLoc;

perFork.doubleForkCount = doubleForkCount;




figure,
histC.minDist = histc(minDist, binning);
plot(binning,histC.minDist/sum(histC.minDist),'-b.')

hold on
histC.minDistImmobilised = histc(minDistImmobilised, binning);
plot(binning,histC.minDistImmobilised/sum(histC.minDistImmobilised),'-r.')


figure,
%histC.minDist = histc(minDist, binning);
plot(binning,histC.minDist,'-b.')

figure,
%histC.minDist = histc(minDist, binning);
plot(binning,histC.minDistImmobilised,'-r.')



figure,
histC.minDistLoc = histc(minDistLoc, binning)';
plot(binning,histC.minDistLoc/sum(histC.minDistLoc),'-b.')

hold on
histC.minDistImmoLoc = histc(minDistImmoLoc, binning)';
plot(binning,histC.minDistImmoLoc/sum(histC.minDistImmoLoc),'-r.')


save([ROIDataFilename '.histC_v2'], 'coLocData');
save([ROIDataFilename '.histCellLength_v2'], 'cellLength')
save([ROIDataFilename '.histPerFork_v2'], 'perFork')

D1c = histc(D1,rangeD'); %normalized histogram count
D2c = histc(D2,rangeD'); %normalized histogram count

figure;
hist(D1,rangeD');
xlim([min(rangeD) max(rangeD)]);
hold all
stem(threshold1,max(D1c),'marker','none');
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');
hold off

figure;
hist(D2,rangeD');
xlim([min(rangeD) max(rangeD)]);
hold all
stem(threshold2,max(D2c),'marker','none');
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');
hold off

diffusionFraction = (ll-1)/(kk-1)

observedFraction = numel(D1)/nMolecules

meanD1 = NaN;
meanD2 = NaN;
countDiscarded
totalForks

%keyboard


end

function [x1 y1] = projectionToLine(imlineVertices,pos)
% Calculate point in predefined line closest to arbitrary point
% Line:     y = m*x + b
% Point:    x0,y0
% Projection in line:   x1,y1

m = diff(imlineVertices(:,2))/diff(imlineVertices(:,1));
b = imlineVertices(1,2)-m*imlineVertices(1,1);
x0 = pos(:,1);
y0 = pos(:,2);

x1 = (m*y0+x0-m*b)/(m^2+1);
y1 = (m^2*y0+m*x0+b)/(m^2+1);
end

