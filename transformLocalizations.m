function transformLocalizations

reply = input('Do you want to generate a tform file? Y/N: ', 's');
if reply == 'Y'
    transformPALMMovies
end

reply = input('Do you want to apply a tform file to localization data? Y/N: ', 's');
if reply == 'Y'
    
    display('Select tform matrix');
    [tformFilename, tformPathname] =...
        uigetfile('*.tform.mat', 'MatFile data:','MultiSelect', 'on');
    
    if ~(isnumeric(tformFilename)&&tformFilename==0) %check the user has not pressed cancel
        
        tformmatrix = importdata([tformPathname tformFilename]);
        
        
        display('Select localization data');
        [dataFilename, dataPathname] =...
            uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
        
        if ~(isnumeric(dataFilename)&&dataFilename==0) %check the user has not pressed cancel
            
            info = whos('dataFilename');
            if strcmp(info.class,'char')
                nFiles = 1;
            else
                nFiles = numel(dataFilename);
            end
            
            for ii = 1:nFiles
                
                if nFiles == 1
                    loadname  = [dataPathname dataFilename];
                else
                    loadname = [dataPathname dataFilename{1,ii}];
                end
                
                newData = importdata(loadname);
                if isstruct(newData)
                    data = newData.data;
                else
                    data = newData;
                end
                clear newData;
                
                
                
                greenPos = data(:,2:3);
                
                bluePos = tforminv(tformmatrix,greenPos);
                
                tformData = data;
                
                tformData(:,2:3) = bluePos;
                
                save([loadname '.tform.out'], 'tformData');
                display(strcat('Saved green - blue transformed localizations as', [loadname '.tform.out']));
                
                clear tformData greenPos bluePos data
                
            end
            
        end
        
    end
    
    
end


end