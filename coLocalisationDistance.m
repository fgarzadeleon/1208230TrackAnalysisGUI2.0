function coLocalisationDistance( pos, ROIData, ROIDataFilename, params )
%COLOCALISATIONDISTANCE Calculates the colocalised distances of tracks and
%position
%   Calculate the distances of the inserted tracks against the position.


set(0,'DefaultLineMarkerSize',5,'DefaultLineLineWidth',2);

cmap = colormap(jet);
pixel = params.pixel;
dT = params.dT;
sigmaNoise = params.sigmaNoise;
threshold = params.DhistThresh;
minSteps = params.DhistMinSteps;
rangeD = params.rangeD;
holdFigure = params.holdFigureCheckbox;
minDist = [];
minDistImmobilised = [];
minDistImmoLoc = [];
minDistLoc = [];

binning = 0:.05:1;
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

for jj = 1:length(ROIData)
    
    impolyVertices = ROIData(jj,1).impolyVertices;
    posINindex = inpolygon(pos(:,1),pos(:,2),impolyVertices(:,1),impolyVertices(:,2));
    posIN = pos(posINindex,:)*pixel;
    
    if sum(posINindex)>2
        doubleForkCount = doubleForkCount+1;
    end
        
    tracksIN = ROIData(jj,1).tracks;
    nMolecules = max(tracksIN(:,4));
    tracksIN(:,1:2) = tracksIN(:,1:2)*pixel;
    plot(pos(:,1), pos(:,2),'ro','MarkerSize',10)
    if any(posINindex~=0)
        
        for ii = 1:nMolecules
            tempDistImmoLoc = [];
            tempDistLoc = [];
            
            xx = find(tracksIN(:,4)==ii);
            
            if numel(xx)>minSteps
                
                MSD = 0; % reset MSD
                
                % sum all squared displacement
                for mm = 1:numel(xx)-1
                    
                    MSD = MSD + ((tracksIN(xx(mm+1),1) - tracksIN(xx(mm),1))^2 +...
                        (tracksIN(xx(mm+1),2) - tracksIN(xx(mm),2))^2);
                    
                end
                
                MSD = MSD/jj; % mean square displacement
                
                % calculate D from MSD and correct for localization noise
                D(kk) = MSD/(4*dT) - sigmaNoise^2*pixel^2/dT;
                
                if D(kk) < threshold
                    
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
                    
                    plot([posIN(indexMin,1) meanPos(ll,1)]/pixel,...
                        [posIN(indexMin,2) meanPos(ll,2)]/pixel,...
                        'b-')
                    
                    plot(meanPos(ll,1)/pixel,meanPos(ll,2)/pixel,'b.')
                    
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
                
                
            end
            kk = kk + 1;
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
histC.minDistLoc = histc(minDistLoc, binning)';
plot(binning,histC.minDistLoc/sum(histC.minDistLoc),'-b.')

hold on
histC.minDistImmoLoc = histc(minDistImmoLoc, binning)';
plot(binning,histC.minDistImmoLoc/sum(histC.minDistImmoLoc),'-r.')


save([ROIDataFilename '.histC'], 'coLocData');
save([ROIDataFilename '.histCellLength'], 'cellLength')
save([ROIDataFilename '.histPerFork'], 'perFork')
end

