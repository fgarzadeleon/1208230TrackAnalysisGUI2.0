function ROIData = selectCellROIs(data,params)


if params.holdFigureCheckbox;
    figure(params.figureHandle);
    hold all;
else
    figure;
end

%initialise data structure
tempStruct = ...
    struct(   'localizationData',zeros(1,2),...
    'tracks',zeros(1,2),...
    'nMolecules',zeros(1,1),...
    'impolyVertices',zeros(1,2),...
    'imlineVertices',zeros(1,2));

remainingData = data;
ROIData = [];

reply = 'Y';
ii = 1;

while reply == 'Y'
    
    v = axis;
    
    h = plot(remainingData(:,2),remainingData(:,3),'b.','MarkerSize',1);
    axis(v);
    
    p = impoly;
    impolyVertices = getPosition(p);
    delete(p);
    
    line = imline;
    input('Done? Hit Return');
    imlineVertices = getPosition(line);
    delete(line);
    
    inIndexes = find( inpolygon(remainingData(:,2),remainingData(:,3),impolyVertices(:,1),impolyVertices(:,2)) );
    inData = remainingData(inIndexes,:);
    
    remainingData(inIndexes,:) = [];
    
    inData(:,11) = ii * ones(numel(inIndexes),1); % ROI id
    
    %track data for this ROI
    if ~isempty(inData)
        pos(:,1) = inData(:,2);
        pos(:,2) = inData(:,3);
        pos(:,3) = inData(:,1);
    else 
        pos = [];
    end
    
    if ~isempty(pos)
        tracks = trackWithDummy(pos, params.trackParams);
        nMolecules = max(tracks(:,4));
        clear pos;
        
        tempStruct.localizationData = inData;
        tempStruct.tracks = tracks;
        tempStruct.nMolecules = nMolecules;
        
    end
    
    tempStruct.impolyVertices = impolyVertices;
    tempStruct.imlineVertices = imlineVertices;
    
    ROIData = [ROIData; tempStruct];
    
    reply = input('Do you want more? Y/N: ', 's');
    if isempty(reply)
        reply = 'Y';
    end
    
    ii = ii + 1;
    delete(h);
    plot([impolyVertices(:,1); impolyVertices(1,1)],[impolyVertices(:,2); impolyVertices(1,2)],'r-.')
    
end

hold off;

%track all data in field of view
pos(:,1) = data(:,2);
pos(:,2) = data(:,3);
pos(:,3) = data(:,1);

tracks = trackWithDummy(pos, params.trackParams);
nMolecules = max(tracks(:,4));
clear pos;

dataStruct.localizationData = data;
dataStruct.tracks = tracks;
dataStruct.nMolecules = nMolecules;
dataStruct.ROIData = ROIData;

end