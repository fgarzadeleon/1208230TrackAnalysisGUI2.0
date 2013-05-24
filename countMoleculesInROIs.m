function countMolecules = countMoleculesInROIs
% load one or many ROI localization files '.ROIData'
% quantify molecules per ROI

% load multiple ROI data
[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'on');

if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
    
    info = whos('ROIDataFilename');
    if strcmp(info.class,'char')
        appData.nFiles = 1;
    else
        appData.nFiles = numel(ROIDataFilename);
    end
    
    kk = 1;
 
    for ii = 1:appData.nFiles
        
        if appData.nFiles == 1
            loadname  = [ROIDataPathname ROIDataFilename];
        else
            loadname = [ROIDataPathname ROIDataFilename{1,ii}];
        end
        
        ROIData = importdata(loadname);
        
        for jj = 1:length(ROIData)
            
            countMolecules(kk) = ROIData(jj,1).nMolecules;
            
            kk = kk + 1;
            
        end
        
    end
            
end

figure;
hist(countMolecules,4:8:150)
v = axis
axis([-20 150 v(3) v(4)+2])
mean(countMolecules)

end