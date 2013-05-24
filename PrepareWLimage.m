function imageI = PrepareWLimage()
%prepare WL images to be read by Schnitzcell (segmoviephase--> FileName-p-001.tif)
%could read .fits and .tif files of same size
%for reading .fits it needs: fits_read_header, ImageStack, getFrame (should work, when gaussStormGUI is working)
%Anne Plochowietz 02/05/12
%implemeted in photobleachingCurveAnalysisGUI
%Anne Plochowietz 11/06/12
%Modified by Federico on 31/08/2012 :)

FileNameCore = 'TempFile';
[WLFilename, WLPathname] =...
    uigetfile({'*.fits';'*.tif'}, 'Brightfield image data:','MultiSelect', 'on');
for fileNum=1:length(WLFilename)
    fileSetName = WLFilename(1:end-5);
    
    if ~(isnumeric(WLFilename)&&WLFilename==0) %check the user has not pressed cancel
        
        info = whos('WLFilename');
        if strcmp(info.class,'char')
            nFiles = 1;
        else
            nFiles = numel(WLFilename);
        end
        
        for ii = 1:nFiles
            
            if nFiles == 1
                loadname  = [WLPathname WLFilename];
            else
                loadname = [WLPathname WLFilename{1,ii}];
            end
            

                ImageInfo = fits_read_header(loadname);
                imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
                
                DataIn=ImageStack(loadname, imageLim);
                imageI = getFrame(DataIn,1);
                %imageI = flipud(imageI); %only needed in connection with
                %tracking

            imageI = uint16(imageI);
            imageI = flipud(imageI);
            
            %image cut depending on boundaries
            
            
            imageI = imcomplement(imageI);
            if (ii<10)  %to load 1 to 99 WL images at once
                FileName = [FileNameCore '-p-00' num2str(ii) '.tif'];
            else
                FileName = [FileNameCore '-p-0' num2str(ii) '.tif'];
            end
            imwrite(imageI,FileName,'tiff');
            imageI = [];    %although loaded images should have the same size (required for Schnitzcell)
        end
    end
end