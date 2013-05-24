function DThresholdInROIs(params)
% load one or many ROI localization files '.ROI.out'
% quantify molecules per ROI

% load and combine multiple ROI data
[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*ROI.out', 'MatFile data:','MultiSelect', 'on');

if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
    
    info = whos('ROIDataFilename');
    if strcmp(info.class,'char')
        appData.nFiles = 1;
    else
        appData.nFiles = numel(ROIDataFilename);
    end
    
    combinedData = [];
    
    for ii = 1:appData.nFiles
        
        if appData.nFiles == 1
            loadname  = [ROIDataPathname ROIDataFilename];
        else
            loadname = [ROIDataPathname ROIDataFilename{1,ii}];
        end
        
        newData = importdata(loadname);
        
        if ii>1
            newData(:,11) = newData(:,11) + max(combinedData(:,11)); %increase the ROI index
        end
        
        combinedData = [combinedData; newData];
        
    end
    
    
else
    
    return
    
end

% track within each ROI
nROIs = max(combinedData(:,11));
countMolecules = zeros(nROIs,1);
diffusionFraction = zeros(nROIs,1);
D1 = zeros(nROIs,1);
D2 = zeros(nROIs,1);

for ii = 1:nROIs
    
    ROIindices = find(combinedData(:,11)==ii);
    
    pos(:,1) = combinedData(ROIindices,2);
    pos(:,2) = combinedData(ROIindices,3);
    pos(:,3) = combinedData(ROIindices,1);
    
    tracks = trackWithDummy(pos, params.trackParams);
    countMolecules(ii) = max(tracks(:,4));
    
    clear pos;
    
    [diffusionFraction(ii), D1(ii), D2(ii)] = twoSpeciesMSDThreshold(tracks, params);
    
end

figure;
hist(diffusionFraction,0:0.1:1)
xlabel('proportion of molecules with D < 0.2');
xlim([0 1]);

figure;
hist(D1);
xlabel('slow diffusion coefficient');

figure;
hist(D2);
xlabel('fast diffusion coefficient');

keyboard;

end