function combineTracksFromLoColi()

combinedTracks = [];

[tracksFilename, tracksPathname] = ...
    uigetfile('*locoli*', 'MatFile data:','MultiSelect', 'on');

if ~(isnumeric(tracksFilename)&&tracksFilename==0) %check the user has not pressed cancel

info = whos('tracksFilename');
if strcmp(info.class,'char')
    appData.nFiles = 1;
else
    appData.nFiles = numel(tracksFilename);
end

for ii = 1:appData.nFiles
    
    if appData.nFiles == 1
        loadname  = [tracksPathname tracksFilename];
    else
        loadname = [tracksPathname tracksFilename{1,ii}];
    end
    
    newTracks = importdata(loadname);
    
    for kk = 1:length(newTracks.cellROI_data)
        ROItracks = newTracks.cellROI_data(kk,1).tracks;
        
        if ii>1 && kk>1
        
            ROItracks(:,4) = ROItracks(:,4) + max(combinedTracks(:,4));
        end
    
    combinedTracks = [combinedTracks; ROItracks];
    end
    
end

[filename, pathname] = uiputfile('.tracks');
save([pathname filename],'combinedTracks');

else
    
    return
    
end

end