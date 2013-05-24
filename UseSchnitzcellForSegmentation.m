function impolyVerticesSchnitzcell = UseSchnitzcellForSegmentation(FileNameCore);
%use Schnitzcell for e.coli WL images from HD and Zugspitze (taken with Andor software)
%date must be in the format 2012-05-02 (2nd of May 2012)
%Anne Plochowietz 02/05/12

%Initialization Schnitzcell
%get current folder
currentfolder = pwd;
%use current folder to save segmentation results

%get current date
date = datestr(now, 'yyyy-mm-dd');


    
    p = initschnitz(FileNameCore,date,'e.coli','rootDir',currentfolder,'imageDir',currentfolder);
    %directory of segmentation results (e.g. FileNameCoreseg001.mat):
    %    currentfolder\date\FileNameCore\segmentation
    
    %Segmentation
    %runs segmentation for all files in format FileNameCore-p-001.tif (increasing numbers)
    %could take up to 1 minute per FOV
    p = segmoviephase(p);
    
    %Manual adjustment of segmentation
    p = manualcheckseg(p);
    
    %get Number of analysed FOVs
    listOfSegmentationFiles = dir([currentfolder '\' date '\' FileNameCore '\segmentation\*.mat']);
    nFiles = numel(listOfSegmentationFiles);
    
    %get segmentation info of single cell per FOV
    for ina1 =1:nFiles %for all FOVs
        %get filenames ready to load
        if(ina1<10)
            %segmentation data
            loadname = [currentfolder '\' date '\' FileNameCore '\segmentation\' FileNameCore 'seg00' num2str(ina1) '.mat'];
            %FL&WL
            %FLname = [FileNameCore '-y-00' num2str(ina1) '.tif'];
            WLname = [FileNameCore '-p-00' num2str(ina1) '.tif'];
        else
            %segmentation data
            loadname = [currentfolder '\' date '\' FileNameCore '\segmentation\' FileNameCore 'seg0' num2str(ina1) '.mat'];
            %FL&WL
            %FLname = [FileNameCore '-y-0' num2str(ina1) '.tif'];
            WLname = [FileNameCore '-p-0' num2str(ina1) '.tif'];
        end
        %load segmentation data into workspace
        segFileData = load(loadname);
        %get number of segmented cells in FOV(ina1)
        numCells = max(max(segFileData.Lc));
        %segFileData.Lc is mask of FOV with numbers 1 to numCells for each single segmented cell
        Cells = segFileData.Lc;
        [r,c] = size(Cells);
        
        %find and load first FLimage, here only first FLimage is considered
%         for ina2=1:1%NumCheckFrames (define number of frames to check for the brightest)
%             Isumsum(1,ina2) = sum(sum(double(imread(FLname,ina2))));
%         end
%         [Isumsummax, firstFrame] = max(Isumsum(1,:));
%         firstFLimage = double(imread(FLname, firstFrame));
%         
        %generate mask of single cell in each field of view with value 1 to get
        %intensity values of the cell
        cellPosition = zeros(numCells,2);
        for ina3 = 1:numCells
            %get mask of single cell
            SingleCell = zeros(r,c);
            SingleCell(Cells == ina3) = 1;
            [indC,indR] = find(SingleCell == 1);
            cellPosition(ina3,:) = [max(indC),max(indR)]; %number of cell will be displayed in the bottom right
            %get intensity values from fluorescence image
            NumPixels = sum(sum(SingleCell));
%             ItotMatrix(ina3,ina1) = sum(sum(firstFLimage(SingleCell==1)));  %total Intensity of Cell
%             IpPMatrix(ina3,ina1) = sum(sum(firstFLimage(SingleCell==1)))/NumPixels;    %Intensity per pixel
        end
        SingleCell = [];
        %open FOV and display Cell numbers or Intensity values on specific cells
        handleFigure = figure;
        set(gcf,'Color',[1 1 1]);
        hold on;
        WLimage = imread(WLname);
        figure(handleFigure), subplot(1,2,1,'position',[0.04 0.1 0.45 0.80]), imagesc(WLimage), colormap(gray), axis equal; hold on;
%         figure(handleFigure), subplot(1,2,2,'position',[0.54 0.1 0.45 0.80]), imagesc(firstFLimage), colormap(gray), axis equal, hold on;
        %cmap = colormap(jet);
        for ina4=1:numCells
            %number of cell and cell itself in FL image
            text(cellPosition(ina4,2),cellPosition(ina4,1),num2str(ina4),'FontSize',14,'Color','y'); hold on;
            %get mask of single cell
            SingleCell = zeros(r,c);
            SingleCell(Cells == ina4) = 1;
            SingleCelledge = edge(SingleCell);
            [y{ina4},x{ina4}] = find(SingleCelledge==1);
            plot(x{ina4},y{ina4},'y.','MarkerSize',6); hold on;
            
            impolyVerticesSchnitzcell{ina4} = [x{ina4},y{ina4}];
            %Color
            %colorindex = ceil( length(cmap) * (ina4) / numCells );
            %plot(...,'Color',cmap(colorindex,:)); hold on;
        end
        hold off;
%         print(['-f' num2str(handleFigure)],'-dtiffn','-r300',[FLname(1:end-4) '_segCells.tif']);
        close(figure(handleFigure));
        
        %draw mask of cells in WL image and in FL image
        h2 = figure;
        set(gcf,'Color',[1 1 1]);
        hold on;
        imshow(WLimage,[min(min(WLimage)) max(max(WLimage))]); hold on;
        for ina5=1:numCells
            plot(x{ina5},y{ina5},'y.','MarkerSize',6); hold on;
        end
        print(['-f' num2str(h2)],'-dtiffn','-r300',[WLname(1:end-4) '_WLonly.tif']);
        close(figure(h2));
        h3 = figure;
        set(gcf,'Color',[1 1 1]);
        hold on;
%         imshow(firstFLimage,[min(min(firstFLimage)) max(max(firstFLimage))]); hold on;
        for ina6=1:numCells
            plot(x{ina6},y{ina6},'y.','MarkerSize',6); hold on;
        end
%         print(['-f' num2str(h3)],'-dtiffn','-r300',[FLname(1:end-4) '_FLonly.tif']);
        close(figure(h3));
        
    end

%store intensity values in textfile
% save([FileNameCore '_IpP.txt'],'IpPMatrix','-ASCII');
% save([FileNameCore '_Itot.txt'],'ItotMatrix','-ASCII');

