function localizationStd = localizationProjection

[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'off');

if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
    
    loadname  = [ROIDataPathname ROIDataFilename];
    
    ROIData = importdata(loadname);
    
    
    kk = 1;
    alignedLocalizations = [];
    
    for jj = 1:length(ROIData)
        
        cellData = ROIData(jj,1).localizationData;
        
        linePositions = ROIData(jj,1).imlineVertices;
        
        x1 = linePositions(1,1);
        y1 = linePositions(1,2);
        x2 = linePositions(2,1);
        y2 = linePositions(2,2);
        
        alpha = atan((y2 - y1)/(x2 - x1));
        
        rotateCellData = cellData;
        rotateCellData(:,2) = cos(alpha)*(cellData(:,2)-x1) + sin(alpha)*(cellData(:,3)-y1);
        rotateCellData(:,3) = -sin(alpha)*(cellData(:,2)-x1) + cos(alpha)*(cellData(:,3)-y1);
        
        % scale cell length to unity
        lineLength(kk) = sqrt( (x2-x1)^2 + (y2-y1)^2 );
        rotateCellData(:,2) = rotateCellData(:,2) / lineLength(kk);
        
        %histogram along line profile
%         figure;
%         hist(rotateCellData(:,2),0:1.6/lineLength(kk):1);
%         xlim([0 1])
        
        localizationMean(kk) = mean(rotateCellData(:,2));
        localizationStd(kk) = std(rotateCellData(:,2));
        
        alignedLocalizations = [alignedLocalizations; rotateCellData(:,2) - localizationMean(kk)];
        
        kk = kk+1;
        
    end
    
    figure;
    hist(alignedLocalizations,-0.5:0.05:0.5);
    xlim([-0.5 0.5])
    xlabel('long cell axis (normalized');
    ylabel('localizations');
    
    figure;
    hist(lineLength);
    xlabel('cell length');
    ylabel('cell count');
    
    figure;
    hist(localizationMean-0.5,-0.5:0.1:0.5);
    xlabel('localization center of mass along long cell axis (normalized)');
    ylabel('cell count');

    
end




% if params.holdFigureCheckbox;
%     f = figure(params.figureHandle);
%     hold all;
% else
%     f = figure;
% end
%
% plot(data(:,2),data(:,3),'.','MarkerSize',3.0);
% xlabel('x [pixels]');
% ylabel('y [pixels]');
%
% reply = 'y';
% kk = 1;
%
% alignedLocalizations = [];
%
% while reply == 'y'
%
%     figure(f);
%     p = impoly;
%     vertices = getPosition(p);
%
%     inIndexes = find( inpolygon(data(:,2),data(:,3),vertices(:,1),vertices(:,2)) );
%     cellData = data(inIndexes,:);
%
%     % draw line for profile
%     line = imline;
%
%     input('Done? Hit Return');
%
%     linePositions = getPosition(line);
%     hold off;
%
%     x1 = linePositions(1,1);
%     y1 = linePositions(1,2);
%     x2 = linePositions(2,1);
%     y2 = linePositions(2,2);
%
%     alpha = atan((y2 - y1)/(x2 - x1));
%
%     rotateCellData = cellData;
%     rotateCellData(:,2) = cos(alpha)*(cellData(:,2)-x1) + sin(alpha)*(cellData(:,3)-y1);
%     rotateCellData(:,3) = -sin(alpha)*(cellData(:,2)-x1) + cos(alpha)*(cellData(:,3)-y1);
%
%     % scale cell length to unity
%     lineLength(kk) = sqrt( (x2-x1)^2 + (y2-y1)^2 );
%     rotateCellData(:,2) = rotateCellData(:,2) / lineLength(kk);
%
%     %histogram along line profile
%     f = figure;
%     hist(rotateCellData(:,2),0:1/lineLength:1);
%     xlim([0 1])
%
%     reply = input('Do you want more? y/n: ', 's');
%     if isempty(reply)
%         reply = 'y';
%     end
%
%     close(f);
%
%     localizationMean(kk) = mean(rotateCellData(:,2));
%     localizationStd(kk) = std(rotateCellData(:,2))
%
%     alignedLocalizations = [alignedLocalizations; rotateCellData(:,2) - localizationMean];
%
%     kk = kk+1;
%
% end
%
% figure;
% hist(alignedLocalizations,-0.5:0.05:0.5);
% xlim([-0.5 0.5])
% xlabel('long cell axis (normalized');
% ylabel('localizations');
%
% figure;
% hist(lineLength);
% xlabel('cell length');
% ylabel('cell count');
%
% figure;
% hist(localizationMean-0.5,-0.5:0.1:0.5);
% xlabel('localization center of mass along long cell axis (normalized)');
% ylabel('cell count');
%
% keyboard;

end
