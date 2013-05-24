function localizationStd = localizationProjection_v2(params)

%Parameters needed
pixel = params.pixel;
minSteps = params.DhistMinSteps;
dT = params.dT;
sigmaNoise = params.sigmaNoise;
threshold1 = params.DThresh(1);
threshold2 = params.DThresh(2);


% Load ROIData
figure
[ROIDataFilename, ROIDataPathname] = ...
    uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'off');

if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
    
    loadname  = [ROIDataPathname ROIDataFilename];
    ROIData = importdata(loadname);
    
    %Load position data
    [posOutFilename, posOutPathname] = ...
        uigetfile('*.out', 'MatFile data:','MultiSelect', 'off');
    
    if ~(isnumeric(posOutFilename)&&posOutFilename==0) %check the user has not pressed cancel
        
        loadname  = [posOutPathname posOutFilename];
        
        posOut = importdata(loadname);
        
        if isfield(posOut,'data')
            
            pos(:,1) = posOut.data(:,2);
            pos(:,2) = posOut.data(:,3);
            allCellData = posOut.data;
            
        else
            allCellData = posOut;
            pos(:,1) = posOut(:,2);
            pos(:,2) = posOut(:,3);
        end
        
        posOut = [];
        
        kk = 1;
        zz = 1;
        alignedLocalizations = [];
        
        locCellData_2 = [];
        locCellData_3 = [];
        
        
        
        for jj = 1:length(ROIData)
            
            posINindex = inpolygon(pos(:,1),pos(:,2),...
                ROIData(jj,1).impolyVertices(:,1),...
                ROIData(jj,1).impolyVertices(:,2));
            
            cellData = [];
            locCellData_2 = [];
            
            locCellData_3 = [];
            
            if 0
                
                tracksIN = ROIData(jj,1).tracks;
                nMolecules = max(tracksIN(:,4));
                tracksIN(:,1:2) = tracksIN(:,1:2)*pixel;
                
                for ii = 1:nMolecules
                    tempDistImmoLoc = [];
                    tempDistLoc = [];
                    
                    xx = find(tracksIN(:,4)==ii);
                    zz = 1;
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
                        D1(zz) = MSD1/(4*dT) - sigmaNoise^2*pixel^2/dT;
                        D2(zz) = MSD2/(8*dT) - sigmaNoise^2*pixel^2/dT;
                        
                        if D1(zz) < threshold1 && D2(zz) < threshold2
                            
                            locCellData_2 = [tracksIN(xx,1)/pixel; locCellData_2];
                            locCellData_3 = [tracksIN(xx,2)/pixel; locCellData_3];
                            
                            %locCellData_2 = [mean(tracksIN(xx,1)/pixel); locCellData_2];
                            %locCellData_3 = [mean(tracksIN(xx,2)/pixel); locCellData_3];
                            
                        end
                        
                        
                        % reset MSD and update index
                        
                        zz = zz + 1;
                    end
                    
                end
                
                
                
                if ~isempty(locCellData_2)
                    
                    cellData(:,2) = locCellData_2;
                    cellData(:,3) = locCellData_3;
                    
                end
            else
                
                cellData = allCellData(posINindex,:);
                
            end
            
            
            if ~isempty(cellData)
                
                linePositions = ROIData(jj,1).imlineVertices;
                
                x1 = linePositions(1,1);
                y1 = linePositions(1,2);
                x2 = linePositions(2,1);
                y2 = linePositions(2,2);
                
                
                
                plot(cellData(:,2), cellData(:,3),'r.')
                hold on
                plot([x1 x2],[y1 y2])
                
                              
                alpha = atan((y2 - y1)/(x2 - x1));
                
                rotateCellData = cellData;
                rotateCellData(:,2) = cos(alpha)*(cellData(:,2)-x1) + sin(alpha)*(cellData(:,3)-y1);
                rotateCellData(:,3) = -sin(alpha)*(cellData(:,2)-x1) + cos(alpha)*(cellData(:,3)-y1);
                
                % scale cell length to unity
                lineLength(kk) = sqrt( (x2-x1)^2 + (y2-y1)^2 );
                rotateCellData(:,2) = rotateCellData(:,2) / lineLength(kk);
                
                if any(rotateCellData(:,2)>2)
                    pause
                end
                
                %histogram along line profile
                %         figure;
                %         hist(rotateCellData(:,2),0:1.6/lineLength(kk):1);
                %         xlim([0 1])
                
                localizationMean(kk) = mean(rotateCellData(:,2));
                localizationStd(kk) = std(rotateCellData(:,2));
                
                cellMeanPos = sqrt(mean([y2-y1,y1-y1])^2+mean([x2-x1,x1-x1])^2);
                
                %alignedLocalizations = [alignedLocalizations; rotateCellData(:,2) - localizationMean(kk)];
                alignedLocalizations = [alignedLocalizations; rotateCellData(:,2)];
                
                
                localizationMean(kk) = localizationMean(kk)-cellMeanPos;
                kk = kk+1;
            end
            
        end
        
        figure;
        hist(alignedLocalizations-0.5,-0.5:.05:0.5);
        xlim([-0.55 0.55])
        xlabel('long cell axis (normalized');
        ylabel('localizations');
        
        figure;
        hist(lineLength);
        xlabel('cell length');
        ylabel('cell count');
        
        figure;
        hist(localizationMean,-0.5:0.1:0.5);
        xlim([-0.5 0.5])
        xlabel('localization center of mass along long cell axis (normalized)');
        ylabel('cell count');
        
        save([posOutFilename '.locProj'], 'alignedLocalizations');
    end
    
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
