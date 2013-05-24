function ROIData = selectROIs(data,params)


if params.holdFigureCheckbox;
    figure(params.figureHandle);
    hold all;
else
    figure;
end

ROIData = [];
remainingData = data;

reply = 'Y';
ii = 1;

while reply == 'Y'
    
    v = axis;

    h = plot(remainingData(:,2),remainingData(:,3),'.');
    axis(v);
    
    p = impoly;
    vertices = getPosition(p);
    delete(p);
    
    inIndexes = find( inpolygon(remainingData(:,2),remainingData(:,3),vertices(:,1),vertices(:,2)) );
    inData = remainingData(inIndexes,:);
    
    remainingData(inIndexes,:) = [];
    
    inData(:,11) = ii * ones(numel(inIndexes),1); % ROI id
    
    ROIData = [ROIData; inData];
    
    reply = input('Do you want more? Y/N: ', 's');
    if isempty(reply)
        reply = 'Y';
    end
    
    ii = ii + 1;
    delete(h);
    
end

hold off;

end
