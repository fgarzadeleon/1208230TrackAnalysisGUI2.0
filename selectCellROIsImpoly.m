function [ROIData, name] = selectCellROIsImpoly(params)


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

%remainingData = data;
ROIData = [];

[name, path] = uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'on');


    ROIDataPre = importdata([path,name]);
    
    for ROIsNum = 1:length(ROIDataPre)
        
        %remainingData = ROIData(ROIsNum).localizationData;
        
        %inIndexes = find( inpolygon(remainingData(:,2),remainingData(:,3),impolyVertices(:,1),impolyVertices(:,2)) );
        %inData = remainingData(inIndexes,:);
        
        %remainingData(inIndexes,:) = [];
        inData = ROIDataPre(ROIsNum).localizationData;
        %inData(:,11) = ii * ones(numel(inIndexes),1); % ROI id
        
        % track data for this ROI
        pos(:,1) = inData(:,2);
        pos(:,2) = inData(:,3);
        pos(:,3) = inData(:,1);
        
        tracks = trackWithDummy(pos, params.trackParams);
        nMolecules = max(tracks(:,4));
        clear pos;
        
        tempStruct.localizationData = inData;
        tempStruct.tracks = tracks;
        tempStruct.nMolecules = nMolecules;
        tempStruct.impolyVertices = ROIDataPre(ROIsNum).impolyVertices;
        tempStruct.imlineVertices = ROIDataPre(ROIsNum).imlineVertices;
        
        ROIData = [ROIData; tempStruct];
        
        %     reply = input('Do you want more? Y/N: ', 's');
        %     if isempty(reply)
        %         reply = 'Y';
        %     end
        
        %     ii = ii + 1;
        %     delete(h);
        
    end

hold off;

% track all data in field of view
% pos(:,1) = data(:,2);
% pos(:,2) = data(:,3);
% pos(:,3) = data(:,1);
%
% tracks = trackWithDummy(pos, params.trackParams);
% nMolecules = max(tracks(:,4));
% clear pos;

% dataStruct.localizationData = data;
% dataStruct.tracks = tracks;
% dataStruct.nMolecules = nMolecules;
% dataStruct.ROIData = ROIData;

end