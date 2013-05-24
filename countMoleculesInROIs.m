function countMolecules = countMoleculesInROIs(appData)
% load one or many ROI localization files '.ROIData'
% quantify molecules per ROI

% load multiple ROI data
[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'on');

nPositionsInCluster = [];
            
distancesWithinClusters = [];

distancesBetweenClusters = [];

if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
    
    info = whos('ROIDataFilename');
    if strcmp(info.class,'char')
        appData.nFiles = 1;
    else
        appData.nFiles = numel(ROIDataFilename);
    end
    
    kk = 1;
 
    for ii = 1:appData.nFiles
        ii
        if appData.nFiles == 1
            loadname  = [ROIDataPathname ROIDataFilename];
        else
            loadname = [ROIDataPathname ROIDataFilename{1,ii}];
        end
        
        ROIData = importdata(loadname);
        
        for jj = 1:length(ROIData)
            
            countMolecules(kk) = ROIData(jj,1).nMolecules;
            
            kk = kk + 1;
            
            pos = zeros(max(ROIData(jj,1).nMolecules),2);
            
            for mols=1:ROIData(jj,1).nMolecules
                
                xx = find(ROIData(jj,1).tracks(:,4)==mols);
                if ~isnan(xx)
                pos(mols,1) = mean(ROIData(jj,1).tracks(xx,1));
                pos(mols,2) = mean(ROIData(jj,1).tracks(xx,2));
                end
                
            end
            
%             clusterThreshold = 5*0.04/appData.pixel;
%             
% 
%             [nPositionsInClusterTemp, distancesWithinClustersTemp,...
%                 distancesBetweenClustersTemp] = nearestNeighbourClustering(pos,...
%                 clusterThreshold, appData);
% 
% 
%             nPositionsInCluster = [nPositionsInCluster; nPositionsInClusterTemp];
% 
%             
%             distancesWithinClusters = [distancesWithinClusters; distancesWithinClustersTemp];
%             
%             distancesBetweenClusters = [distancesBetweenClusters; distancesBetweenClustersTemp];
            
        end
        
    end
            
end

figure; hist(nPositionsInCluster,20);
xlabel('Numbers')
title('Number of positions within clusters');

figure; hist(distancesWithinClusters,100);
xlabel('distance')
title('Distances within clusters');

figure; hist(distancesBetweenClusters,10);
xlabel('distance')
title('Distances between clusters');


figure;
hist(countMolecules,4:8:150)
v = axis
axis([-20 150 v(3) v(4)+2])
mean(countMolecules)
median(countMolecules)


end